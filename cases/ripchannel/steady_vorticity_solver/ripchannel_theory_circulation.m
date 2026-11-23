function out = ripchannel_theory_circulation(info,opts)
%RIPCHANNEL_THEORY_CIRCULATION Wave-averaged theoretical rip-current maps.
%
% out = ripchannel_theory_circulation(info,opts)
%
% INPUT
%   info.x, info.y : coordinates [m], x cross-shore (shore -> offshore),
%                    y alongshore. Vectors.
%   info.h         : still-water depth [m], size [numel(y),numel(x)],
%                    positive in wet cells. This is the convention used by
%                    make_ripchannel_bathy.m.
%   info.Hs        : offshore significant wave height [m].
%   info.Tp        : peak period [s].
%   info.spread    : directional standard deviation [degrees].
%   info.Dp        : mean direction [degrees] relative to shore-normal.
%
% OPTIONS (all optional)
%   opts.source       'bonneton' (default), 'bowen', or 'both'.
%   opts.theta0       override info.Dp [degrees]. Positive points toward +y.
%   opts.spread       override info.spread [degrees].
%   opts.ndir         number of directional bins; default  nine (odd).
%   opts.gamma        depth-limited breaker ratio Hb/h; default 0.78.
%   opts.nu           horizontal eddy viscosity [m^2/s]; default 5.
%   opts.drag         linear vorticity damping rate [1/s]; default 0.03.
%   opts.hmin         minimum wet depth [m]; default 0.25.
%   opts.maxCells     maximum analysis grid cells; default 180000.
%   opts.breakWidth   smooth breaking transition, as a fraction of Hb;
%                    default 0.08. Set to 0 for a hard threshold.
%   opts.doPlot       make diagnostic figures; default true.
%   opts.saveFile     MAT-file path; empty means do not save.
%
% OUTPUT
%   out.x,out.y,out.h                  analysis grid and depth
%   out.Hrms,out.Dissipation           directional-integrated wave fields
%   out.vorticitySource                Bonneton curl(D e_k) [1/s^2]
%   out.vorticitySourceBowen           curl of breaking-zone radiation-stress
%                                      acceleration [1/s^2]
%   out.omega,out.psi                  modelled relative vorticity and
%                                      transport streamfunction
%   out.Ux,out.Uy                      depth-averaged currents [m/s]
%   out.directional                    directional components and weights
%
% THEORY / IMPORTANT SCOPE
%   Bonneton et al. (2010), equations (8)-(9) and (18),(22), motivate the
%   dissipative force D and vorticity source curl(D e_k). Equation (8) is
%   used here as a saturated-breaking parameterization, component by
%   component. The shock expression (18) requires bore states h1,h2 and
%   bore celerity cb; those are not available from info.h alone, so this
%   implementation does NOT claim to evaluate (18) exactly.
%
%   Bowen (1969) motivates a radiation-stress forcing confined to the surf
%   zone and a steady circulation response to longshore breaker-height
%   variations. Here, the radiation-stress divergence is evaluated from a
%   directional spectrum and its curl is used as the vorticity forcing.
%   The steady response is obtained from a linear vorticity diffusion-drag
%   closure and a Poisson inversion for the transport streamfunction. This
%   is a practical 2-D numerical extension, not Bowen's analytic sinusoidal
%   solution or a complete nonlinear shallow-water model.
%
%   Directional spreading is represented by independent direction bins
%   with energy weights, linear shoaling and Snell refraction over the local
%   cross-shore depth. This first implementation does not solve a full
%   2-D spectral-action/ray equation, so alongshore refraction/focusing of
%   each direction is not explicitly resolved. Use it as a theoretical
%   baseline to compare with FUNWAVE-TVD, not as a replacement for it.
%
% Coordinates: x increases offshore. Incident wave propagation is toward
% shore, so e_k = [-cos(theta), sin(theta)]. This direction convention is
% stated explicitly because FUNWAVE's direction convention should be
% checked against the values in info.Dp for a given case.

if nargin < 2 || isempty(opts), opts = struct(); end
opts = setdef(opts,'source','bonneton');
opts = setdef(opts,'theta0',getfielddef(info,'Dp',0));
opts = setdef(opts,'spread',getfielddef(info,'spread',0));
opts = setdef(opts,'ndir',9);
opts = setdef(opts,'gamma',0.78);
opts = setdef(opts,'nu',5.0);
opts = setdef(opts,'drag',0.03);
opts = setdef(opts,'hmin',0.25);
opts = setdef(opts,'maxCells',180000);
opts = setdef(opts,'breakWidth',0.08);
opts = setdef(opts,'doPlot',true);
opts = setdef(opts,'saveFile','');

required = {'x','y','h','Hs','Tp'};
for n = 1:numel(required)
    if ~isfield(info,required{n})
        error('info.%s is required. Load it using prep_local_ripchannel_info and make_ripchannel_bathy.',required{n});
    end
end
x0 = info.x(:)'; y0 = info.y(:); h0 = info.h;
if ~isequal(size(h0),[numel(y0),numel(x0)])
    error('info.h must have size [numel(info.y), numel(info.x)].');
end
if numel(x0)<3 || numel(y0)<3, error('At least 3 points are required in x and y.'); end
if any(diff(x0)<=0) || any(diff(y0)<=0), error('info.x and info.y must be strictly increasing.'); end
if info.Hs<=0 || info.Tp<=0, error('info.Hs and info.Tp must be positive.'); end
if opts.nu<=0 || opts.drag<=0, error('opts.nu and opts.drag must be positive.'); end

% Coarsen large FUNWAVE grids for a tractable sparse elliptic solve.
stride = max(1,ceil(sqrt(numel(h0)/opts.maxCells)));
ix = 1:stride:numel(x0);
iy = 1:stride:numel(y0);
x = x0(ix); y = y0(iy); h = double(h0(iy,ix));
[XX,YY] = meshgrid(x,y); %#ok<ASGLU>
wet = isfinite(h) & h>opts.hmin;
h(~wet) = opts.hmin;

% Directional quadrature: Gaussian in direction, normalized to unit energy.
sigma = max(0,double(opts.spread));
if sigma < 1e-8
    theta = deg2rad(opts.theta0);
    weights = 1;
else
    nd = max(3,round(opts.ndir)); if mod(nd,2)==0, nd=nd+1; end
    thmax = min(3*sigma,75);
    thetaDeg = linspace(-thmax,thmax,nd);
    thetaDeg = thetaDeg + opts.theta0;
    theta = deg2rad(thetaDeg);
    weights = exp(-0.5*((thetaDeg-opts.theta0)/sigma).^2);
    weights = weights/sum(weights);
end

% Linear dispersion and group speed on the analysis grid.
g = 9.81; rho = 1025; omega = 2*pi/info.Tp;
hOff = median(h(:,end));
hOff = max(hOff,opts.hmin);
[kOff,cOff,cgOff] = dispersion(omega,hOff,g);

ndir = numel(theta); ny = numel(y); nx = numel(x);
Dsum = zeros(ny,nx); FxD = zeros(ny,nx); FyD = zeros(ny,nx);
Sxx = zeros(ny,nx); Sxy = zeros(ny,nx); Syy = zeros(ny,nx);
H2sum = zeros(ny,nx); breakingEnergy = zeros(ny,nx);
Hcomp = zeros(ny,nx,ndir); Dcomp = zeros(ny,nx,ndir);
for q = 1:ndir
    th0 = theta(q); wq = weights(q);
    % Each bin carries its fraction of wave variance.
    H0q = info.Hs*sqrt(wq);
    [k,c,cg] = dispersion(omega,h,g);
    % Snell refraction over cross-shore depth contours.
    sinTh = sin(th0).*c./cOff;
    sinTh = max(-0.995,min(0.995,sinTh));
    th = asin(sinTh);
    cth = max(cos(th),0.05);
    cth0 = max(cos(th0),0.05);
    % Conservation of shore-normal wave-energy flux before breaking.
    Hraw = H0q.*sqrt((cgOff*cth0)./(max(cg.*cth,1e-6)));
    Hlim = opts.gamma.*h;
    ratio = Hraw./max(Hlim,1e-6);
    if opts.breakWidth>0
        br = 0.5*(1+tanh((ratio-1)/opts.breakWidth));
    else
        br = double(ratio>=1);
    end
    Hloc = min(Hraw,Hlim);
    % Bonneton et al. (2010), equation (8), component-wise saturated form.
    Dq = br.*g.*Hloc.^3 ./ (4*max(c,0.05).*info.Tp.*max(h,opts.hmin).^2);
    % Wave propagation vector (x is positive offshore).
    ex = -cos(th); ey = sin(th);
    Dsum = Dsum + Dq;
    FxD = FxD + Dq.*ex;
    FyD = FyD + Dq.*ey;
    H2sum = H2sum + Hloc.^2;
    breakingEnergy = breakingEnergy + br.*(Hloc.^2);
    Hcomp(:,:,q) = Hraw;
    Dcomp(:,:,q) = Dq;

    % Directional radiation-stress tensor, Longuet-Higgins form.
    E = rho*g*Hloc.^2/8;
    n = cg./max(c,0.05);
    Sxx = Sxx + E.*(n.*ex.^2 + (n-0.5));
    Sxy = Sxy + E.*(n.*ex.*ey);
    Syy = Syy + E.*(n.*ey.^2 + (n-0.5));
end
Hrms = sqrt(H2sum);

% Bonneton et al. (2010), equation (9)/(22): curl(D e_k).
[dFy_dx,~] = gradient(FyD,dx,dy);
[~,dFx_dy] = gradient(FxD,dx,dy);
sourceBonneton = dFy_dx - dFx_dy;
sourceBonneton(~wet) = 0;

% Bowen-type radiation-stress forcing. Set-up is irrotational in the
% idealized plane-beach/offshore region; the surf-zone stress divergence
% is retained as the forcing for this reduced vorticity calculation.
[dSxx_dx,dSxx_dy] = gradient(Sxx,dx,dy);
[dSxy_dx,dSxy_dy] = gradient(Sxy,dx,dy);
FxRad = -(dSxx_dx + dSxy_dy)./(rho.*max(h,opts.hmin));
% Correct y-component: dSxy/dx + dSyy/dy.
[~,dSyy_dy] = gradient(Syy,dx,dy);
FyRad = -(dSxy_dx + dSyy_dy)./(rho.*max(h,opts.hmin));
[~,dFxRad_dy] = gradient(FxRad,dx,dy);
[dFyRad_dx,~] = gradient(FyRad,dx,dy);
sourceBowen = dFyRad_dx - dFxRad_dy;
% Retain forcing only where at least one direction is breaking.
sourceBowen = sourceBowen.*double(breakingEnergy>0);
sourceBowen(~wet) = 0;

% Select the source. 'both' returns two independent model solutions.
sourceName = lower(char(opts.source));
switch sourceName
    case 'bonneton'
        [omegaMap,psi,Ux,Uy] = solve_circulation(sourceBonneton,h,x,y,opts);
    case 'bowen'
        [omegaMap,psi,Ux,Uy] = solve_circulation(sourceBowen,h,x,y,opts);
    case 'both'
        [omegaMap,psi,Ux,Uy] = solve_circulation(sourceBonneton,h,x,y,opts);
        [omegaBowen,psiBowen,UxBowen,UyBowen] = solve_circulation(sourceBowen,h,x,y,opts);
    otherwise
        error('opts.source must be ''bonneton'', ''bowen'', or ''both''.');
end

out = struct();
out.x=x; out.y=y; out.h=h; out.wet=wet;
out.Hrms=Hrms; out.Dissipation=Dsum;
out.Fd_x=FxD; out.Fd_y=FyD;
out.vorticitySource=sourceBonneton;
out.vorticitySourceBowen=sourceBowen;
out.Sxx=Sxx; out.Sxy=Sxy; out.Syy=Syy;
out.omega=omegaMap; out.psi=psi; out.Ux=Ux; out.Uy=Uy;
out.directional.theta=theta; out.directional.weights=weights;
out.directional.H=Hcomp; out.directional.D=Dcomp;
out.options=opts; out.coarsenStride=stride;
if strcmp(sourceName,'both')
    out.bowen.omega=omegaBowen; out.bowen.psi=psiBowen;
    out.bowen.Ux=UxBowen; out.bowen.Uy=UyBowen;
end

if opts.doPlot
    plot_results(out,sourceName);
end
if ~isempty(opts.saveFile)
    save(opts.saveFile,'out','-v7.3');
end
end

function [om,psi,Ux,Uy] = solve_circulation(src,h,x,y,opts)
% Steady linear vorticity diffusion-drag closure:
%   0 = nu Laplacian(omega) - r omega + source
% followed by -Laplacian(psi)=omega and M=(dpsi/dy,-dpsi/dx).
% r is a linear damping rate (Bowen's linear bottom-friction analogue).
ny=numel(y); nx=numel(x); N=ny*nx;
dx=median(diff(x)); dy=median(diff(y));
ex=ones(nx,1); ey=ones(ny,1);
Tx=spdiags([ex -2*ex ex],[-1 0 1],nx,nx)/dx^2;
Ty=spdiags([ey -2*ey ey],[-1 0 1],ny,ny)/dy^2;
L=kron(speye(nx),Ty)+kron(Tx,speye(ny));
% Fixed zero-vorticity and zero-streamfunction perimeter. This makes a
% closed theoretical domain; enlarge the domain to reduce boundary effects.
bd=false(ny,nx); bd([1 end],:)=true; bd(:,[1 end])=true; bd=bd(:);
R=opts.drag*ones(ny,nx); R(~isfinite(h))=opts.drag;
A= -opts.nu*L + spdiags(R(:),0,N,N);
A=dirichlet(A,bd);
b=src(:); b(bd)=0;
omv=solve_spd(A,b); om=reshape(omv,ny,nx);
Apsi=dirichlet(-L,bd);
psiv=solve_spd(Apsi,om(:)); psi=reshape(psiv,ny,nx);
[dpsi_dx,dpsi_dy]=gradient(psi,dx,dy);
Mx=dpsi_dy; My=-dpsi_dx;
Ux=Mx./max(h,opts.hmin); Uy=My./max(h,opts.hmin);
% Mask dry cells and a narrow numerical boundary rim.
wet=h>opts.hmin; Ux(~wet)=NaN; Uy(~wet)=NaN; om(~wet)=NaN; psi(~wet)=NaN;
end

function A=dirichlet(A,bd)
idx=find(bd); A(idx,:)=0; A(:,idx)=0; A(idx,idx)=speye(numel(idx));
end

function x=solve_spd(A,b)
if norm(b)==0, x=zeros(size(b)); return; end
try
    setup.type='ict'; setup.droptol=2e-3; setup.diagcomp=1e-3;
    L=ichol(A,setup);
    [x,flag]=pcg(A,b,1e-7,1200,L,L');
catch
    flag=1; x=zeros(size(b));
end
if flag~=0
    % A sparse direct solve is a robust fallback for modest grids.
    if size(A,1)<=250000
        x=A\b;
    else
        warning('PCG did not fully converge (flag %d). Returning its iterative result.',flag);
    end
end
end

function [k,c,cg]=dispersion(w,h,g)
% Vectorized Newton iteration for w^2=g*k*tanh(k*h).
h=max(h,0.05); k=max(w^2/g,w./sqrt(g*h));
for it=1:30
    kh=k.*h; th=tanh(kh); sech2=1-th.^2;
    f=g*k.*th-w^2;
    fp=g*(th+k.*h.*sech2);
    dk=f./max(fp,1e-12);
    k=max(1e-8,k-dk);
    if max(abs(dk(:)))<1e-10, break; end
end
c=w./k; kh=k.*h; z=2*kh;
fac=ones(size(z)); ii=z>1e-6; fac(ii)=z(ii)./sinh(min(z(ii),700));
fac(z>700)=0;
cg=0.5*c.*(1+fac);
end

function plot_results(out,sourceName)
[xg,yg]=meshgrid(out.x,out.y);
figure('Color','w','Name','Wave-averaged rip-current theory');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
contourf(xg,yg,-out.h,24,'LineStyle','none'); axis ij tight; colorbar;
xlabel('Cross-shore x [m]'); ylabel('Alongshore y [m]'); title('Bed elevation z_b=-h [m]');
nexttile;
if strcmp(sourceName,'bowen'), q=out.vorticitySourceBowen; else, q=out.vorticitySource; end
contourf(xg,yg,q,25,'LineStyle','none'); axis ij tight; colorbar;
xlabel('Cross-shore x [m]'); ylabel('Alongshore y [m]'); title('Wave-induced vorticity source [s^{-2}]');
nexttile;
pcolor(xg,yg,sqrt(out.Ux.^2+out.Uy.^2)); shading flat; hold on;
step=max(1,round(min(numel(out.x),numel(out.y))/35));
quiver(xg(1:step:end,1:step:end),yg(1:step:end,1:step:end), ...
    out.Ux(1:step:end,1:step:end),out.Uy(1:step:end,1:step:end),1.5,'k');
contour(xg,yg,out.psi,12,'k','LineWidth',0.65);
axis ij tight; colorbar; xlabel('Cross-shore x [m]'); ylabel('Alongshore y [m]');
title('Current speed, vectors and transport streamlines');
sgtitle(sprintf('Wave-averaged circulation: %s source; spread = %.1f deg',sourceName,out.options.spread));
if strcmp(sourceName,'both')
    figure('Color','w','Name','Bowen-type circulation');
    pcolor(xg,yg,sqrt(out.bowen.Ux.^2+out.bowen.Uy.^2)); shading flat; hold on;
    step=max(1,round(min(numel(out.x),numel(out.y))/35));
    quiver(xg(1:step:end,1:step:end),yg(1:step:end,1:step:end), ...
        out.bowen.Ux(1:step:end,1:step:end),out.bowen.Uy(1:step:end,1:step:end),1.5,'k');
    contour(xg,yg,out.bowen.psi,12,'k','LineWidth',0.65);
    axis ij tight; colorbar; xlabel('Cross-shore x [m]'); ylabel('Alongshore y [m]');
    title('Bowen-type radiation-stress current speed and streamlines');
end
end

function s=setdef(s,name,value)
if ~isfield(s,name) || isempty(s.(name)), s.(name)=value; end
end
function v=getfielddef(s,name,default)
if isfield(s,name) && ~isempty(s.(name)), v=s.(name); else, v=default; end
end
