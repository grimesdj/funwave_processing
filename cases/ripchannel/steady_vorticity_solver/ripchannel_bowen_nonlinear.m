function out = ripchannel_bowen_nonlinear(info,opts)
%RIPCHANNEL_BOWEN_NONLINEAR Nonlinear, wave-averaged rip-current model.
%
% out = ripchannel_bowen_nonlinear(info,opts)
%
% Solves the steady depth-averaged shallow-water equations in pseudo-time:
%
%   d(eta)/dtau + div(h U) = 0
%   dU/dtau + U.grad(U) = -g grad(eta)
%       - (1/(rho h)) div(S) - cL U/h + nu laplacian(U),
%
% with h = d + eta. The implementation is advanced in conservative
% perturbation variables using a Rusanov (local Lax-Friedrichs) flux. A
% hydrostatic reference-flux subtraction makes the scheme well-balanced for
% a resting state over the prescribed FUNWAVE bathymetry.
%
% This is a numerical nonlinear extension of Bowen (1969), not a verbatim
% reproduction of every approximation in that paper. In particular, it
% uses the directional wave/radiation-stress fields produced by
% ripchannel_theory_circulation.m and solves the steady nonlinear shallow-
% water equations by pseudo-time relaxation. Alongshore boundaries are
% periodic; cross-shore boundaries are closed, with an offshore damping
% sponge to approximate decay away from the surf zone.
%
% REQUIRED INPUT
%   info.x, info.y, info.h, info.Hs, info.Tp
%   info.spread and info.Dp are used when available.
%
% OPTIONS
%   opts.maxCells       maximum working-grid cells (default 80000)
%   opts.gamma          breaker index passed to wave model (default 0.78)
%   opts.ndir           number of directional bins (default 9)
%   opts.friction       'bottom' (default), 'eddy', or 'both'. Bowen
%                       treats linear bottom drag and eddy viscosity as
%                       alternative closures; 'both' is a sensitivity option.
%   opts.nu             horizontal eddy viscosity [m^2/s], default 2
%   opts.bottomDrag     Bowen linear drag parameter cL [m/s], default 0.02
%   opts.hmin           minimum wet depth [m], default 0.20
%   opts.cfl             pseudo-time CFL number, default 0.30
%   opts.maxIter         maximum iterations, default 60000
%   opts.checkEvery      convergence check interval, default 100
%   opts.tol             normalized steady residual tolerance, default 1e-6
%   opts.relax           update under-relaxation, default 0.8
%   opts.spongeFraction  offshore sponge fraction of x-domain, default 0.12
%   opts.spongeStrength  maximum sponge rate [1/s], default 0.05
%   opts.doPlot          make plots, default true
%   opts.saveFile        MAT-file output path, default ''
%   opts.verbose         print iteration diagnostics, default true
%
% OUTPUT
%   out.x,out.y,out.d,out.h,out.eta,out.Ux,out.Uy,out.omega
%   out.Sxx,out.Sxy,out.Syy,out.Fwave_x,out.Fwave_y
%   out.iterations,out.converged,out.residualHistory,out.options
%
% The stress fields are calculated by ripchannel_theory_circulation.m.
% Keep that function on the MATLAB path in the same directory as this file.
%
% Coordinates follow the uploaded make_ripchannel_bathy.m convention:
% x increases offshore; y is periodic alongshore. Depth d is positive.
%
% References:
%   Bowen, A.J. (1969), Rip currents. 1. Theoretical investigations,
%   J. Geophys. Res., 74, 5467-5478, especially equations (19)-(28).
%   Bonneton, P. et al. (2010), Large-scale vorticity generation due to
%   dissipating waves in the surf zone, DCDS-B, 13, 729-738.

if nargin < 2 || isempty(opts), opts = struct(); end
opts = setdef(opts,'maxCells',80000);
opts = setdef(opts,'gamma',0.78);
opts = setdef(opts,'ndir',9);
opts = setdef(opts,'friction','bottom');
opts = setdef(opts,'nu',2.0);
opts = setdef(opts,'bottomDrag',0.02);
opts = setdef(opts,'hmin',0.20);
opts = setdef(opts,'cfl',0.30);
opts = setdef(opts,'maxIter',60000);
opts = setdef(opts,'checkEvery',100);
opts = setdef(opts,'tol',1e-6);
opts = setdef(opts,'relax',0.8);
opts = setdef(opts,'spongeFraction',0.12);
opts = setdef(opts,'spongeStrength',0.05);
opts = setdef(opts,'doPlot',true);
opts = setdef(opts,'saveFile','');
opts = setdef(opts,'verbose',true);

req = {'x','y','h','Hs','Tp'};
for k=1:numel(req)
    if ~isfield(info,req{k}), error('info.%s is required.',req{k}); end
end
if opts.nu < 0 || opts.bottomDrag < 0 || opts.hmin <= 0
    error('nu and bottomDrag must be nonnegative; hmin must be positive.');
end
friction = lower(char(opts.friction));
switch friction
    case 'bottom', nuEff = 0; dragEff = opts.bottomDrag;
    case 'eddy',   nuEff = opts.nu; dragEff = 0;
    case 'both',   nuEff = opts.nu; dragEff = opts.bottomDrag;
    otherwise, error('opts.friction must be ''bottom'', ''eddy'', or ''both''.');
end
if strcmp(friction,'bottom') && dragEff<=0
    error('bottomDrag must be positive for friction=''bottom''.');
end
if strcmp(friction,'eddy') && nuEff<=0
    error('nu must be positive for friction=''eddy''.');
end

% Generate the same directional wave field and grid coarsening used by the
% preceding linear theory routine. The linear currents are only an initial
% guess; the final solution is from the nonlinear equations below.
wopt.source = 'bowen';
wopt.theta0 = fielddef(info,'Dp',0);
wopt.spread = fielddef(info,'spread',0);
wopt.ndir = opts.ndir;
wopt.gamma = opts.gamma;
wopt.nu = max(opts.nu,0.1);
wopt.drag = 0.03;
wopt.hmin = opts.hmin;
wopt.maxCells = opts.maxCells;
wopt.doPlot = false;
wopt.saveFile = '';
wave = ripchannel_theory_circulation(info,wopt);

x = wave.x(:)'; y = wave.y(:); d = double(wave.h);
[ny,nx] = size(d);
if nx < 4 || ny < 4
    error('Working grid is too small. Increase info grid size or maxCells.');
end
dx = median(diff(x)); dy = median(diff(y));
if abs(max(diff(x))-min(diff(x))) > 1e-6*dx || ...
   abs(max(diff(y))-min(diff(y))) > 1e-6*dy
    error('This solver requires uniformly spaced x and y coordinates.');
end

% Fixed land mask. Wet/dry motion is not included in this steady solver.
wet = isfinite(d) & d > opts.hmin;
d(~wet) = 0;
if nnz(wet) < 20, error('Too few wet cells after applying hmin.'); end

% Bowen radiation-stress force per unit mass. S has units N/m, so
% -div(S)/(rho*h) is a horizontal acceleration.
g = 9.81; rho = 1025;
[Sxx,Sxy,Syy] = deal(wave.Sxx,wave.Sxy,wave.Syy);
Sxx(~wet)=0; Sxy(~wet)=0; Syy(~wet)=0;
divSx = ddx(Sxx,dx) + ddy_periodic(Sxy,dy);
divSy = ddx(Sxy,dx) + ddy_periodic(Syy,dy);
Fx = -divSx ./ (rho .* max(d,opts.hmin));
Fy = -divSy ./ (rho .* max(d,opts.hmin));
Fx(~wet)=0; Fy(~wet)=0;

% Surface elevation and depth-integrated transports Q=(eta,hu,hv).
eta = zeros(ny,nx);
Ux = wave.Ux; Uy = wave.Uy;
Ux(~isfinite(Ux))=0; Uy(~isfinite(Uy))=0;
Ux(~wet)=0; Uy(~wet)=0;
% No-normal flow at shore/offshore ends; the offshore sponge handles the
% transition from the active surf zone to the far-field boundary.
coast = wet & ~circshift(wet,[0 1]);
offshore = false(ny,nx); offshore(:,end)=true;
Ux(coast | offshore)=0; Uy(coast | offshore)=0;

% Sponge ramps smoothly to zero at the offshore end.
ns = max(2,round(opts.spongeFraction*nx));
sponge = zeros(1,nx);
if ns >= 2
    s = linspace(0,1,ns);
    sponge(nx-ns+1:nx) = opts.spongeStrength .* (0.5-0.5*cos(pi*s));
end
sponge = repmat(sponge,ny,1);
sponge(~wet)=0;

% Iteration diagnostics. Residual is the max norm of the pseudo-time
% tendency, normalized by its initial value. A second criterion checks
% relative changes in the free surface and transport.
resHist = nan(ceil(opts.maxIter/opts.checkEvery)+1,1);
dtHist = nan(size(resHist));
res0 = NaN; ncheck=0; converged=false; lastChange=Inf;

% A CFL step based on the maximum shallow-water wave speed and velocity.
% The explicit viscosity restriction is also enforced.
for it = 1:opts.maxIter
    h = zeros(ny,nx);
    h(wet) = max(d(wet)+eta(wet),opts.hmin);
    u = zeros(ny,nx); v = zeros(ny,nx);
    u(wet)=Ux(wet); v(wet)=Uy(wet);
    qx = h.*u; qy = h.*v;

    cmax = max(sqrt(g*h(wet)));
    umax = max(sqrt(u(wet).^2+v(wet).^2));
    dtWave = opts.cfl*min(dx,dy)/max(cmax+umax,1e-6);
    if nuEff > 0
        dtVisc = 0.20/(nuEff*(1/dx^2+1/dy^2));
    else
        dtVisc = Inf;
    end
    dt = min(dtWave,dtVisc);

    % Rusanov fluxes for perturbation variables about the resting state
    % eta=0, U=0 over the actual bathymetry d(x,y). This avoids numerical
    % diffusion of the nonuniform still-water depth.
    [F1,F2,F3,G1,G2,G3] = fluxes_rusanov(eta,qx,qy,h,d,wet,u,v,g,dx,dy);
    divQ = (F1(:,2:end)-F1(:,1:end-1))/dx + ...
           (G1([2:end 1],:)-G1)/dy;
    divMx = (F2(:,2:end)-F2(:,1:end-1))/dx + ...
            (G2([2:end 1],:)-G2)/dy;
    divMy = (F3(:,2:end)-F3(:,1:end-1))/dx + ...
            (G3([2:end 1],:)-G3)/dy;

    ddxBed = ddx(d,dx); ddyBed = ddy_periodic(d,dy);
    lapU = laplacian_periodic_y(u,dx,dy);
    lapV = laplacian_periodic_y(v,dx,dy);
    % Bowen's linear bottom-friction parameterization is proportional to
    % velocity divided by total depth in the primitive momentum equation.
    % Multiplication by h gives -cL*U in the transport equations.
    rhsEta = -divQ;
    rhsMx = -divMx + g.*eta.*ddxBed - divSx./rho ...
            - dragEff.*u + nuEff.*h.*lapU;
    rhsMy = -divMy + g.*eta.*ddyBed - divSy./rho ...
            - dragEff.*v + nuEff.*h.*lapV;

    % Offshore sponge: relaxation toward still water and zero current.
    rhsEta = rhsEta - sponge.*eta;
    rhsMx = rhsMx - sponge.*qx;
    rhsMy = rhsMy - sponge.*qy;

    rhsEta(~wet)=0; rhsMx(~wet)=0; rhsMy(~wet)=0;
    rhsEta(:,[1 end])=0; rhsMx(:,[1 end])=0; rhsMy(:,[1 end])=0;

    % Convert conservative momentum tendencies to velocity tendencies.
    % qx_t = h*u_t + u*h_t, so the quotient rule is required.
    etaDot = rhsEta;
    uDot = (rhsMx - u.*rhsEta)./max(h,opts.hmin);
    vDot = (rhsMy - v.*rhsEta)./max(h,opts.hmin);
    uDot(~wet)=0; vDot(~wet)=0;

    % Residual uses the physical tendencies before under-relaxation.
    rnow = max([max(abs(etaDot(wet))), ...
                max(abs(uDot(wet)))*min(dx,dy)/max(cmax,1e-6), ...
                max(abs(vDot(wet)))*min(dx,dy)/max(cmax,1e-6)]);
    if ~isfinite(rnow), error('Non-finite residual at iteration %d.',it); end

    % Forward pseudo-time step with under-relaxation.
    etaNew = eta + opts.relax*dt*etaDot;
    uNew = u + opts.relax*dt*uDot;
    vNew = v + opts.relax*dt*vDot;
    etaNew(~wet)=0; uNew(~wet)=0; vNew(~wet)=0;
    % Do not allow the fixed wet mask to become shallower than hmin.
    etaNew(wet)=max(etaNew(wet),opts.hmin-d(wet)+1e-5);
    % Landward shoreline is impermeable; impose Bowen's zero-velocity
    % shoreline condition on both components at the first wet cells.
    uNew(coast)=0; vNew(coast)=0;
    % Closed offshore edge and smooth sponge just landward of it.
    uNew(:,end)=0; vNew(:,end)=0; etaNew(:,end)=0;

    lastChange = max([max(abs(etaNew(wet)-eta(wet))), ...
                      max(abs(uNew(wet)-u(wet))), ...
                      max(abs(vNew(wet)-v(wet)))]);
    eta=etaNew; Ux=uNew; Uy=vNew;

    if mod(it,opts.checkEvery)==0 || it==1
        ncheck=ncheck+1;
        if isnan(res0), res0=max(rnow,1e-14); end
        resHist(ncheck)=rnow/res0;
        dtHist(ncheck)=dt;
        if opts.verbose
            fprintf('Bowen nonlinear: iter %6d, residual %.3e, change %.3e, dt %.3e s\n', ...
                it,resHist(ncheck),lastChange,dt);
        end
        if resHist(ncheck)<opts.tol && lastChange<opts.tol*max(1,max(abs([Ux(wet);Uy(wet)])))
            converged=true;
            break
        end
    end
end

% Final depth, wet mask and vorticity.
h = zeros(ny,nx); h(wet)=max(d(wet)+eta(wet),opts.hmin);
Ux(~wet)=NaN; Uy(~wet)=NaN; eta(~wet)=NaN;
uc=Ux; vc=Uy; uc(~wet)=0; vc(~wet)=0;
omega = ddx(vc,dx)-ddy_periodic(uc,dy);
omega(~wet)=NaN;

out = struct();
out.x=x; out.y=y; out.d=d; out.h=h; out.eta=eta; out.wet=wet;
out.Ux=Ux; out.Uy=Uy; out.omega=omega;
out.Sxx=Sxx; out.Sxy=Sxy; out.Syy=Syy;
out.Fwave_x=Fx; out.Fwave_y=Fy;
out.iterations=it; out.converged=converged;
out.residualHistory=resHist(1:ncheck);
out.dtHistory=dtHist(1:ncheck);
out.lastChange=lastChange; out.options=opts;
out.effectiveFriction=friction; out.effectiveNu=nuEff; out.effectiveBottomDrag=dragEff;
out.wave=wave;
if ~converged
    warning('ripchannel_bowen_nonlinear:NoConvergence', ...
        'Solver stopped after %d iterations without reaching tol=%.2e. Inspect residualHistory.',it,opts.tol);
end

if opts.doPlot, plot_nonlinear(out); end
if ~isempty(opts.saveFile), save(opts.saveFile,'out','-v7.3'); end
end

function [F1,F2,F3,G1,G2,G3] = fluxes_rusanov(eta,qx,qy,h,d,wet,u,v,g,dx,dy)
% Rusanov numerical fluxes. x faces are represented by ny-by-(nx+1)
% arrays; y faces are periodic and represented by ny-by-nx arrays.
[ny,nx]=size(h);
% x-face states: left cell i, right cell i+1.
HL=h(:,1:end-1); HR=h(:,2:end);
EL=eta(:,1:end-1); ER=eta(:,2:end);
UL=u(:,1:end-1); UR=u(:,2:end);
VL=v(:,1:end-1); VR=v(:,2:end);
DL=d(:,1:end-1); DR=d(:,2:end);
WL=wet(:,1:end-1); WR=wet(:,2:end);
% Close faces adjacent to land, and at the two cross-shore boundaries.
a=max(abs(UL)+sqrt(g*HL),abs(UR)+sqrt(g*HR));
Q1L=EL; Q1R=ER; Q2L=HL.*UL; Q2R=HR.*UR; Q3L=HL.*VL; Q3R=HR.*VR;
fx1L=HL.*UL; fx1R=HR.*UR;
fx2L=HL.*UL.^2+0.5*g*(HL.^2-DL.^2);
fx2R=HR.*UR.^2+0.5*g*(HR.^2-DR.^2);
fx3L=HL.*UL.*VL; fx3R=HR.*UR.*VR;
fx1=0.5*(fx1L+fx1R)-0.5*a.*(Q1R-Q1L);
fx2=0.5*(fx2L+fx2R)-0.5*a.*(Q2R-Q2L);
fx3=0.5*(fx3L+fx3R)-0.5*a.*(Q3R-Q3L);
faceWet=WL & WR;
fx1(~faceWet)=0; fx2(~faceWet)=0; fx3(~faceWet)=0;
F1=zeros(ny,nx+1); F2=F1; F3=F1;
F1(:,2:nx)=fx1; F2(:,2:nx)=fx2; F3(:,2:nx)=fx3;

% y-face: lower row j, upper row j+1, with periodic wrap.
HU=h; UU=u; VV=v; EE=eta; DD=d;
HR=circshift(HU,[-1 0]); UR=circshift(UU,[-1 0]);
VR=circshift(VV,[-1 0]); ER=circshift(EE,[-1 0]);
DR=circshift(DD,[-1 0]); WR=circshift(wet,[-1 0]);
WL=wet; UL=UU; VL=VV; EL=EE; DL=DD; HL=HU;
a=max(abs(VL)+sqrt(g*HL),abs(VR)+sqrt(g*HR));
Q1L=EL; Q1R=ER; Q2L=HL.*UL; Q2R=HR.*UR; Q3L=HL.*VL; Q3R=HR.*VR;
gy1L=HL.*VL; gy1R=HR.*VR;
gy2L=HL.*UL.*VL; gy2R=HR.*UR.*VR;
gy3L=HL.*VL.^2+0.5*g*(HL.^2-DL.^2);
gy3R=HR.*VR.^2+0.5*g*(HR.^2-DR.^2);
gy1=0.5*(gy1L+gy1R)-0.5*a.*(Q1R-Q1L);
gy2=0.5*(gy2L+gy2R)-0.5*a.*(Q2R-Q2L);
gy3=0.5*(gy3L+gy3R)-0.5*a.*(Q3R-Q3L);
faceWet=WL & WR;
gy1(~faceWet)=0; gy2(~faceWet)=0; gy3(~faceWet)=0;
G1=gy1; G2=gy2; G3=gy3;
end

function z=ddx(a,dx)
% Centered derivative with zero-gradient edge closure.
z=zeros(size(a));
z(:,2:end-1)=(a(:,3:end)-a(:,1:end-2))/(2*dx);
z(:,1)=(a(:,2)-a(:,1))/dx;
z(:,end)=(a(:,end)-a(:,end-1))/dx;
end
function z=ddy_periodic(a,dy)
z=(circshift(a,[-1 0])-circshift(a,[1 0]))/(2*dy);
end
function z=laplacian_periodic_y(a,dx,dy)
z=zeros(size(a));
z(:,2:end-1)=(a(:,3:end)-2*a(:,2:end-1)+a(:,1:end-2))/dx^2;
z(:,1)=(a(:,2)-a(:,1))/dx^2;
z(:,end)=(a(:,end-1)-a(:,end))/dx^2;
z=z+(circshift(a,[-1 0])-2*a+circshift(a,[1 0]))/dy^2;
end
function s=setdef(s,n,v)
if ~isfield(s,n) || isempty(s.(n)), s.(n)=v; end
end
function v=fielddef(s,n,d)
if isfield(s,n) && ~isempty(s.(n)), v=s.(n); else, v=d; end
end
function plot_nonlinear(out)
[XX,YY]=meshgrid(out.x,out.y);
spd=sqrt(out.Ux.^2+out.Uy.^2);
figure('Color','w','Name','Nonlinear Bowen circulation');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
contourf(XX,YY,-out.d,24,'LineStyle','none'); hold on;
contour(XX,YY,out.wet, [0.5 0.5],'k');
axis ij tight; colorbar; xlabel('x [m]'); ylabel('y [m]');
title('Bed elevation z_b=-d [m]');
nexttile;
contourf(XX,YY,out.eta,24,'LineStyle','none');
axis ij tight; colorbar; xlabel('x [m]'); ylabel('y [m]');
title('Mean water-level anomaly [m]');
nexttile;
pcolor(XX,YY,spd); shading flat; hold on;
step=max(1,round(min(numel(out.x),numel(out.y))/35));
quiver(XX(1:step:end,1:step:end),YY(1:step:end,1:step:end), ...
    out.Ux(1:step:end,1:step:end),out.Uy(1:step:end,1:step:end),1.4,'k');
contour(XX,YY,out.eta,12,'k','LineWidth',0.6);
axis ij tight; colorbar; xlabel('x [m]'); ylabel('y [m]');
title('Nonlinear current speed and vectors');
sgtitle(sprintf('Nonlinear Bowen solver: %d iterations, converged = %d', ...
    out.iterations,out.converged));
figure('Color','w','Name','Nonlinear solver convergence');
semilogy(1:numel(out.residualHistory),out.residualHistory,'LineWidth',1.5);
grid on; xlabel('Convergence check'); ylabel('Normalized residual');
title('Pseudo-time convergence history');
end
