function info = estimate_FUNWAVE_rotational_velocity_alongshore_spectra(info);
%
% get the input (x,y,h) grid file
load(info.bathyFile)
% switch from h=-z_bottom to h=depth
% h = -h;
% preserve the original grid
x0=x;
y0=y;
h0=h;
%
% full alongshore domain
if ~isfield(info,'subDomain')
    iX = find(x0>=0  & x0<=400);
    nx = length(iX);
    ny = length(y0);
    subDomain = [1 ny iX(1) iX(end)];
else
    subDomain = info.subDomain;
end
%
% map bottom points to eta points
x = x(subDomain(3):subDomain(4));
y = y(subDomain(1):subDomain(2));
h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
             h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
Nx = length(x);
Ny = length(y);
dy = median(diff(y));
dx = median(diff(x));
%
% i) get the x-coordinates in each bin
ROT = matfile(info.rotVelFile,'writable',true);
t       = ROT.t; t = t-t(1);
Nt      = length(t);
dt      = mean(diff(t));
% for mapping to the local shoreline coordinates:
xsl        = info.x_shoreline';
xp         = x-mean(xsl);
x_map      = x-xsl;
[xx,yy]    = meshgrid(x-mean(xsl),y);
% ii) load the u_rot/v_rot in each bin
Urot     = ROT.Urot;
Vrot     = ROT.Vrot;
VORT     = ROT.VORT;
%
% compute statistics separately for Urot = <Urot> + Urot'
Urot_avg = mean(Urot,3,'omitnan');
Vrot_avg = mean(Vrot,3,'omitnan');
VORT_avg = mean(VORT,3,'omitnan');
%
% 4b) plot y-avg[], y-rms[], and y-spectra[] of Prot_avg
% first, map to shoreline coordinates
if std(xsl)>=3
    fprintf('\tmapping to shoreline coordinates\n')
    Urot_avg_map = griddata(x_map,yy,Urot_avg,xx,yy);
    Vrot_avg_map = griddata(x_map,yy,Vrot_avg,xx,yy);
    VORT_avg_map = griddata(x_map,yy,VORT_avg,xx,yy);
else
    Urot_avg_map = Urot_avg;
    Vrot_avg_map = Vrot_avg;
    VORT_avg_map = VORT_avg;
end
    % estimate spectra
[Urot_avg_spec,ky] = alongshore_spectra_estimate(info,Urot_avg_map);
[Vrot_avg_spec,ky] = alongshore_spectra_estimate(info,Vrot_avg_map);
[VORT_avg_spec,ky] = alongshore_spectra_estimate(info,VORT_avg_map);
Nk = length(ky);
%
Urot_avg_yavg = mean(Urot_avg_map,1,'omitnan');
Vrot_avg_yavg = mean(Vrot_avg_map,1,'omitnan');
VORT_avg_yavg = mean(VORT_avg_map,1,'omitnan');
Urot_yavg     = nan(1,Nx,Nt);
Vrot_yavg     = nan(1,Nx,Nt);
VORT_yavg     = nan(1,Nx,Nt);
Urot_spec     = nan(Nk,Nx,Nt);
Vrot_spec     = nan(Nk,Nx,Nt);
VORT_spec     = nan(Nk,Nx,Nt);
KEflux        = nan(Nx,Nt);
%
fprintf('\tlooping over time to estimage spectra\n')
for jj = 1:Nt
    % 4b) plot y-avg[], y-rms[], and y-spectra[] of Prot_avg
    % first, map to shoreline coordinates if std(xsl)>3m (FRF shorelines can be crazy)
    if std(xsl)>=3
        fprintf('\tmapping to shoreline coordinates\n')
        Urot_map = griddata(x_map,yy,Urot(:,:,jj)-Urot_avg,xx,yy);
        Vrot_map = griddata(x_map,yy,Vrot(:,:,jj)-Vrot_avg,xx,yy);
        VORT_map = griddata(x_map,yy,VORT(:,:,jj)-VORT_avg,xx,yy);
    else
        Urot_map = Urot(:,:,jj)-Urot_avg;
        Vrot_map = Vrot(:,:,jj)-Vrot_avg;
        VORT_map = VORT(:,:,jj)-VORT_avg;
    end
    % estimate spectra
    [Utmp_spec,ky] = alongshore_spectra_estimate(info,Urot_map);
    [Vtmp_spec,ky] = alongshore_spectra_estimate(info,Vrot_map);
    [VVVV_spec,ky] = alongshore_spectra_estimate(info,VORT_map);
    Urot_spec(:,:,jj) = Utmp_spec;
    Vrot_spec(:,:,jj) = Vtmp_spec;
    VORT_spec(:,:,jj) = VVVV_spec;    
    Urot_yavg(1,:,jj) = mean(Urot_map,1,'omitnan');
    Vrot_yavg(1,:,jj) = mean(Vrot_map,1,'omitnan');
    VORT_yavg(1,:,jj) = mean(VORT_map,1,'omitnan');
    KEflux(:,Nt)      = mean((Urot_map+Urot_avg_map).*((Urot_map+Urot_avg_map).^2 + (Vrot_map+Vrot_avg_map).^2),1,'omitnan');
end
%
Urot_yavg = mean(Urot_yavg,3,'omitnan');
Vrot_yavg = mean(Vrot_yavg,3,'omitnan');
VORT_yavg = mean(VORT_yavg,3,'omitnan');
Urot_spec = mean(Urot_spec,3,'omitnan');
Vrot_spec = mean(Vrot_spec,3,'omitnan');
VORT_spec = mean(VORT_spec,3,'omitnan');
ROT.KEflux = mean(KEflux,3,'omitnan');
%
% now cross-shore average into 25m bins (dof~2*50, less the cross-shore correlation)
db    = 50;
xbins = [0:db:150]+db/2;
Nb    = length(xbins);
Urot_spec_binned = [];
Vrot_spec_binned = [];
VORT_spec_binned = [];
Urot_avg_spec_binned = [];
Vrot_avg_spec_binned = [];
VORT_avg_spec_binned = [];
fprintf('\taveraging spectra into cross-shore bins\n')
for ii = 1:Nb
    inbin  =  find( xp>=(xbins(ii)-db/2) & xp<(xbins(ii)+db/2));
    Urot_spec_binned(:,ii) = mean(Urot_spec(:,inbin),2,'omitnan');
    Vrot_spec_binned(:,ii) = mean(Vrot_spec(:,inbin),2,'omitnan');
    VORT_spec_binned(:,ii) = mean(VORT_spec(:,inbin),2,'omitnan');    
    Urot_avg_spec_binned(:,ii) = mean(Urot_avg_spec(:,inbin),2,'omitnan');
    Vrot_avg_spec_binned(:,ii) = mean(Vrot_avg_spec(:,inbin),2,'omitnan');
    VORT_avg_spec_binned(:,ii) = mean(VORT_avg_spec(:,inbin),2,'omitnan');    
end
ROT.ky = ky;
ROT.Urot_yavg = Urot_yavg;
ROT.Vrot_yavg = Vrot_yavg;
ROT.VORT_yavg = VORT_yavg;
ROT.Urot_spec = Urot_spec;
ROT.Vrot_spec = Vrot_spec;
ROT.VORT_spec = VORT_spec;
%
ROT.Urot_spec_binned = Urot_spec_binned;
ROT.Vrot_spec_binned = Vrot_spec_binned;
ROT.VORT_spec_binned = VORT_spec_binned;
%
ROT.Urot_avg_yavg = Urot_avg_yavg;
ROT.Vrot_avg_yavg = Vrot_avg_yavg;
ROT.VORT_avg_yavg = VORT_avg_yavg;
ROT.Urot_avg_spec = Urot_avg_spec;
ROT.Vrot_avg_spec = Vrot_avg_spec;
ROT.VORT_avg_spec = VORT_avg_spec;
%
ROT.Urot_avg_spec_binned = Urot_avg_spec_binned;
ROT.Vrot_avg_spec_binned = Vrot_avg_spec_binned;
ROT.VORT_avg_spec_binned = VORT_avg_spec_binned;
%
% 4a) make 3-panel map of Prot_avg and components,
xm = 2; ym = 2; pw = 10; ph = 3; ag=0.25; ps = [1.5*xm+pw, 2*ym+4*ag+3*ph];
ppos1 = [xm ym               pw ph];
ppos2 = [xm ym+ph+ag         pw ph];
ppos3 = [xm ym+2*ph+2*ag     pw ph];
cbpos = [xm+0.75*pw ym+1.74*ph+2*ag 0.25*pw ag/2];
fig0  = figure('units','centimeters');
pos   = get(fig0,'position');
pos(3:4) = ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
%
% y-avg & y-rms fields
f0ax0  = axes('units','centimeters','position',ppos3);
% $$$ p0 = loglog(f0ax0,ky,mean(Urot_avg_spec_binned,2,'omitnan'),'-r','linewidth',2);
f0ax0.XAxis.Scale='log';
f0ax0.YAxis.Scale='log';
hold(f0ax0,'on')
cm2   = cmocean('thermal',Nb);
for ii = 1:Nb
    loglog(f0ax0,ky,Urot_avg_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
p0 = loglog(f0ax0,ky,mean(Urot_spec_binned,2,'omitnan'),'--r','linewidth',1);
%
ylabel(f0ax0,'$S_{\langle u\rangle,\langle u\rangle}$~[(m/s)$^2$/$\Delta k$]','interpreter','latex','fontsize',12)
% f0l1  = legend(p0,'$S_{u,u}$');
% set(f0l1,'interpreter','latex','fontsize',12)
title(info.runName)
set(f0ax0,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xticklabel',[],'xminortick','on','box','on','ylim',[min(min(Urot_avg_spec_binned(:)),min(Vrot_avg_spec_binned(:)))  max(max(Urot_avg_spec_binned(:)),max(Vrot_avg_spec_binned(:)))])
grid(f0ax0,'on')
%
% y-avg & y-rms fields
f0ax1  = axes('units','centimeters','position',ppos2);
% p1 = loglog(f0ax1,ky,mean(Vrot_avg_spec_binned,2,'omitnan'),'-b','linewidth',2);
f0ax1.XAxis.Scale='log';
f0ax1.YAxis.Scale='log';
hold(f0ax1,'on')
for ii = 1:Nb
    loglog(f0ax1,ky,Vrot_avg_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
p0 = loglog(f0ax1,ky,mean(Vrot_spec_binned,2,'omitnan'),'--r','linewidth',1);
%xlabel(f0ax1,' $k_y$~[1/m]','interpreter','latex')
ylabel(f0ax1,'$S_{\langle v\rangle,\langle v\rangle}$','interpreter','latex','fontsize',12)
%f0l1  = legend(p1,'');
% set(f0l1,'interpreter','latex','fontsize',12)
set(f0ax1,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xticklabel',[],'xminortick','on','box','on','ylim',[min(min(Urot_avg_spec_binned(:)),min(Vrot_avg_spec_binned(:)))  max(max(Urot_avg_spec_binned(:)),max(Vrot_avg_spec_binned(:)))])
grid(f0ax1,'on')
%
f0ax2  = axes('units','centimeters','position',ppos1);
% loglog(f0ax2,ky,mean(VORT_avg_spec_binned,2,'omitnan'),'-k','linewidth',2)
f0ax2.XAxis.Scale='log';
f0ax2.YAxis.Scale='log';
hold(f0ax2,'on')
for ii = 1:Nb
    loglog(f0ax2,ky,VORT_avg_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
%
xlabel(f0ax2,' $k_y$~[1/m]','interpreter','latex')
ylabel(f0ax2,'$S_{\langle \omega\rangle\langle \omega\rangle}$~[(1/s)$^2$/$\Delta k$]','interpreter','latex','fontsize',12)
% f0l2  = legend('$S_{\omega,\omega}$');
% set(f0l2,'interpreter','latex','fontsize',12)
set(f0ax2,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xminortick','on','box','on')
grid(f0ax2,'on')
set(f0ax2,'tickdir','out','ticklabelinterpreter','latex','box','on')
%
f0cb = axes('units','centimeters','position',cbpos);
imagesc(f0cb,xbins,0,reshape(cm2,1,Nb,3))
xlabel(f0cb,'$(x-\bar{x}_\mathrm{sl})$','interpreter','latex','rotation',0,'fontsize',12)
set(f0cb, 'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','bottom','ytick',[],'ticklength',2*get(f0cb,'ticklength'),'xtick',xbins)
%
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_velocity_and_vorticity_time_avg_spectra.pdf'])
%
%
%
%
fig1  = figure('units','centimeters');
pos   = get(fig1,'position');
pos(3:4) = ps;
set(fig1,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
%
% y-avg & y-rms fields
f0ax0  = axes('units','centimeters','position',ppos3);
f0ax0.XAxis.Scale='log';
f0ax0.YAxis.Scale='log';
hold(f0ax0,'on')
cm2   = cmocean('thermal',Nb);
for ii = 1:Nb
    loglog(f0ax0,ky,Urot_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
p0 = loglog(f0ax0,ky,mean(Urot_avg_spec_binned,2,'omitnan'),'--r','linewidth',1);
%
ylabel(f0ax0,'$S_{u,u}$~[(m/s)$^2$/$\Delta k$]','interpreter','latex','fontsize',12)
% f0l1  = legend(p0,'$S_{u,u}$');
% set(f0l1,'interpreter','latex','fontsize',12)
title(info.runName)
set(f0ax0,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xminortick','on','box','on','ylim',[min(min(Urot_spec_binned(:)),min(Vrot_spec_binned(:)))  max(max(Urot_spec_binned(:)),max(Vrot_spec_binned(:)))])
grid(f0ax0,'on')
%
% y-avg & y-rms fields
f0ax1  = axes('units','centimeters','position',ppos2);
f0ax1.XAxis.Scale='log';
f0ax1.YAxis.Scale='log';
hold(f0ax1,'on')
for ii = 1:Nb
    loglog(f0ax1,ky,Vrot_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
p1 = loglog(f0ax1,ky,mean(Vrot_avg_spec_binned,2,'omitnan'),'--r','linewidth',1);
%xlabel(f0ax1,' $k_y$~[1/m]','interpreter','latex')
ylabel(f0ax1,'$S_{v,v}$','interpreter','latex','fontsize',12)
%f0l1  = legend(p1,'');
% set(f0l1,'interpreter','latex','fontsize',12)
set(f0ax1,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xminortick','on','box','on','ylim',[min(min(Urot_spec_binned(:)),min(Vrot_spec_binned(:)))  max(max(Urot_spec_binned(:)),max(Vrot_spec_binned(:)))])
grid(f0ax1,'on')
%
f0ax2  = axes('units','centimeters','position',ppos1);
% loglog(f0ax2,ky,mean(VORT_spec_binned,2,'omitnan'),'-k','linewidth',2)
f0ax2.XAxis.Scale='log';
f0ax2.YAxis.Scale='log';
hold(f0ax2,'on')
for ii = 1:Nb
    loglog(f0ax2,ky,VORT_spec_binned(:,ii),'-','color',cm2(ii,:),'linewidth',2),
end
p1 = loglog(f0ax2,ky,mean(VORT_avg_spec_binned,2,'omitnan'),'--r','linewidth',1);
%
xlabel(f0ax2,' $k_y$~[1/m]','interpreter','latex')
ylabel(f0ax2,'$S_{\omega\omega}$~[(1/s)$^2$/$\Delta k$]','interpreter','latex','fontsize',12)
% f0l2  = legend('$S_{\omega,\omega}$');
% set(f0l2,'interpreter','latex','fontsize',12)
set(f0ax2,'tickdir','out','ticklabelinterpreter','latex','xlim',[1e-3 1],'xminortick','on','box','on')
grid(f0ax2,'on')
set(f0ax2,'tickdir','out','ticklabelinterpreter','latex','box','on')
%
f0cb = axes('units','centimeters','position',cbpos);
imagesc(f0cb,xbins,0,reshape(cm2,1,Nb,3))
xlabel(f0cb,'$(x-\bar{x}_\mathrm{sl})$','interpreter','latex','rotation',0,'fontsize',12)
set(f0cb, 'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','bottom','ytick',[],'ticklength',2*get(f0cb,'ticklength'),'xtick',xbins)
%
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_velocity_and_vorticity_spectra.pdf'])
%
%
%
stationsFile = [info.rootMat,info.rootName,'gauges.mat'];
if exist(stationsFile,'file')
    load(stationsFile)
    stationsX = cell2mat({stations.x}');
    stationsY = cell2mat({stations.y}');
    stationsUrot = cell2mat({stations.Urot}');
else
    stationsX = [];
    stationsY = [];
    stationsUrot = [];
end
%
%
fig2  = figure('units','centimeters');
pos   = get(fig2,'position');
pos(3:4) = ps;
set(fig2,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
% y-avg & y-rms fields
f2ax0  = axes('units','centimeters','position',ppos3);
dk = ky(2)-ky(1);
p0 = plot(f2ax0,xp,sqrt(sum(Urot_avg_spec,1)*dk),'-r',xp,sqrt(sum(Vrot_avg_spec,1)*dk),'-b',stationsX,stationsUrot,'sk','linewidth',2);
%
ylabel(f2ax0,'$\sigma$~[m/s]','interpreter','latex','fontsize',12)
f2l1  = legend(p0(1:2),'$\langle u\rangle$','$\langle v\rangle$');
set(f2l1,'interpreter','latex','fontsize',12)
title(info.runName)
set(f2ax0,'tickdir','out','ticklabelinterpreter','latex','xlim',-25+[0 300],'xticklabel',[])
grid(f2ax0,'on')
%
% y-avg & y-rms fields
f2ax1  = axes('units','centimeters','position',ppos2);
p1 = plot(f2ax1,xp,sqrt(sum(Urot_spec,1)*dk),'-r',xp,sqrt(sum(Vrot_spec,1)*dk),'-b',stationsX,stationsUrot,'sk','linewidth',2);
ylabel(f2ax1,'$\sigma$~[m/s]','interpreter','latex','fontsize',12)
f2l1  = legend(p0(1:2),'$u''$','$v''$');
set(f2l1,'interpreter','latex','fontsize',12)
set(f2ax1,'tickdir','out','ticklabelinterpreter','latex','xlim',-25+[0 300],'xticklabel',[])
grid(f2ax1,'on')
%
f2ax2  = axes('units','centimeters','position',ppos1);
plot(f2ax2,xp,sqrt(sum(VORT_spec,1)*dk),'-k','linewidth',2)
hold(f2ax2,'on')
xlabel(f2ax2,' $x$~[m]','interpreter','latex')
ylabel(f2ax2,'$\sigma_{\omega}$~[(1/s)]','interpreter','latex','fontsize',12)
% f2l2  = legend('$S_{\omega,\omega}$');
% set(f2l2,'interpreter','latex','fontsize',12)
set(f2ax2,'tickdir','out','ticklabelinterpreter','latex','xlim',-25+[0 300])
grid(f2ax2,'on')
set(f2ax2,'tickdir','out','ticklabelinterpreter','latex','box','on')
%
exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_velocity_and_vorticity_standard_deviation.pdf'])
%
