function info = estimate_FUNWAVE_spectral_energy_flux(info);
%
% USAGE: info = estimate_FUNWAVE_spectral_energy_flux(info);
%
% estimate the spectral energy flux in Fourier/wavenumber space,
% beginning with the spectral transport term,
%           T(k)   =  ~(u_rot)* ~((u_rot\dot\nabla) u_rot),
% where u_rot is the rotational term of the helmholtz decompostion,
%           u      = u_rot + u_irr.
% The flux is estimated based on the integrated transport due to all
% wavenumbers >=k,
%           Pi(k)  = cumsum( dk * T(k), 'reverse' ).
% Analysis is limited to the cross-shore region:
%       shoreline <= x <=  525 m

% $$$ % get the input (x,y,h) grid file
% $$$ load(info.bathyFile)
% load grid(h,x,y) and subDomain used to make u_rot
load(info.waveStatsFile,'h','x','y','subDomain')
% switch from h=-z_bottom to h=depth
h = -h;
% preserve the original grid
x0=x;
y0=y;
h0=h;
%
% $$$ % map bottom points to eta points
% $$$ x = x(subDomain(3):subDomain(4));
% $$$ y = y(subDomain(1):subDomain(2));
% $$$ h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
% $$$ h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
% $$$              h(2:end,1:end-1) + h(2:end,2:end));
% $$$ x  = 0.5*(x(1:end-1) + x(2:end));
% $$$ y  = 0.5*(y(1:end-1) + y(2:end));
%
nx = length(x);
ny = length(y);
dy = median(diff(y));
dx = median(diff(x));
%
% map to a shoreline coordinate system
xsl       = info.x_shoreline;
xsl_avg   = mean(xsl);
x_map     = x-xsl';
%
% get pointer to data-file
ROT    = matfile(info.rotVelFile);
t      = ROT.t;
nt     = length(t);
Ns     = 300;
Ne     = floor(nt/Ns);
%
% for averages
inds   = 0;
T      = 0;
P      = 0;
%
for ii = 1:(Ne-1)
    % define indices for this burst
    inds = inds(end) + [1:Ns];
    % load velocity
    u    = ROT.Urot(:,:,inds);
    v    = ROT.Vrot(:,:,inds);
    %
    % estimate y-advection term: (u v_x + v v_y)
    dvdx = (v(:,3:end,:)-v(:,1:end-2,:))/(2*dx);
    dvdy = (v(3:end,:,:)-v(1:end-2,:,:))/(2*dy);    
    % make gradient=0 at boundaries
    dvdx = cat(2,dvdx(:,1,:), dvdx, dvdx(:,end,:));
    dvdy = cat(1,(v(2,:,:)-v(end,:,:))/(2*dy), dvdy, (v(1,:,:)-v(end-1,:,:))/(2*dy));    
    %
    ay_x = u.*dvdx;
    ay_y = v.*dvdy;
    %
    % Fourier transform all fields
    [Ay_x,ky] = alongshore_fft_estimate(info,ay_x);
    [Ay_y, ~] = alongshore_fft_estimate(info,ay_y);
    [V   , ~] = alongshore_fft_estimate(info,v);
    dk        = ky(2)-ky(1);
    %    [U   , ~] = alongshore_fft_estimate(info,u);
    %
    Tx = -conj(V).*Ay_x;
    Ty = -conj(V).*Ay_y;
    T  = T + mean(Tx+Ty,3,'omitnan');    
    %
    Px = cumsum(Tx*dk,1,'reverse');
    Py = cumsum(Ty*dk,1,'reverse');
    P  = P + mean(Px+Py,3,'omitnan');
end
%
T = T/Ne;
P = P/Ne;
%
% cross-shore bin average
dx_bin = 25;
x_bins = [floor(xsl_avg)+dx_bin]:dx_bin:max(x)-dx_bin;
N_bins = length(x_bins);
Pavg   = [];
Tavg   = [];
for jj = 1:length(x_bins)
    inds = find(x>=x_bins(jj)-dx_bin/2 & x<x_bins(jj)+dx_bin/2);
    Pavg = cat(2,Pavg, mean(P(:,inds),2,'omitnan'));
    Tavg = cat(2,Tavg, mean(T(:,inds),2,'omitnan'));    
end
%
%
% plot
xm = 2;
ym = 2;
pw = 8;
ph = 5;
ag = 0.25;
ppos = [xm ym pw ph];
cbpos= [xm+pw+ag, ym, ag, ph/2];
ps   = [2*xm+pw+2*ag, 2*ym+ph];
%
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
cm = cmocean('thermal',N_bins);
ax  = axes('units','centimeters','position',ppos);
semilogx(ax,ky,real(Pavg))
grid on
xlabel(ax,'$k_y$ [m$^{-1}$]','interpreter','latex')
ylabel(ax,'$\Pi(k_y,x)$ [m$^{2}$s$^{-3}/k_y$]','interpreter','latex')
title(info.runName)
set(ax,'ColorOrder',cm,'xtick',[1e-3 1e-2 1e-1],'xlim',[1e-3 1e-1],...
       'tickdir','out','ticklabelinterpreter','latex')
%
% get title string
% $$$ subdirs = split(info.rootSim,'/');
% $$$ runID   = subdirs{end-1};
% $$$ title(ax,runID,'interpreter','latex')
%
%
cb  = axes('units','centimeters','position',cbpos);
imagesc(cb,0,(x_bins-xsl_avg)/xsl_avg,reshape(cm,[N_bins,1,3]))
xlabel(cb,{'~~~$x/L_\mathrm{\scriptscriptstyle{SZ}}$';'~~~~[m]'},'interpreter','latex')
set(cb,'tickdir','out','xtick',[],'xaxislocation','top','yaxislocation','right',...
       'ticklabelinterpreter','latex','ydir','normal')
%
figname = [info.rootSim,filesep,'figures',filesep,info.rootName,'spectral_energy_flux.png'];
exportgraphics(fig,figname)
%
close all
info.specFluxFile = [info.rootMat,info.rootName,'spectral_energy_flux.mat'];
save(info.specFluxFile,'ky','x','P','T','x_bins','Pavg','Tavg')
%
save(info.fileName,'-struct','info')
