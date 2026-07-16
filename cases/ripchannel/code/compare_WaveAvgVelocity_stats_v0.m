%% create compendium plot of exchange and enstrophy for a set of runs:
clear all
close all
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'

%% size for height
possible_heights = [0.5 1.0 1.5];
marker_sizes     = [6 8 10];

%% colormap for directional spread:
possible_spreads = [0 2 4 10 20];
N = length(possible_spreads);
cm = cmocean('thermal',N+1);
cm = cm(1:N,:);

%% markers for bathymetry:
possible_bathy = {'barRip0','barRip1','terRip1'};
bathy_markers  = {'o','s','d'};

%% 1) need a list of run-directories
NAMES    = {'highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10','spreadRip-barRip0','spreadRip-barRip1','spreadRip-terRip1'};

fig0 = figure;
fig1 = figure;
fig2 = figure;
for nn = 1:length(NAMES)
    
    NAME     = NAMES{nn}
switch NAME
  case 'uniRip-ter2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'};
    run_dirs = {'ter2D_h10t10s02d00','ter2D_h10t10s04d00','ter2D_h10t10s10d00','ter2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'uniRip-bar2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'};
    run_dirs = {'bar2D_h10t10s02d00','bar2D_h10t10s04d00','bar2D_h10t10s10d00','bar2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'highRip-barRip0-s00'
    runIDs   = {'highRip','spreadRip','highRip'};
    run_dirs = {'barRip0_h05t10s00d00','barRip0_h10t10s00d00','barRip0_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'};   
    run_dirs = {'barRip1_h05t10s00d00','barRip1_h10t10s00d00','barRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'};   
    run_dirs = {'terRip1_h05t10s00d00','terRip1_h10t10s00d00','terRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip0-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'barRip0_h05t10s10d00','barRip0_h10t10s10d00','barRip0_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'barRip1_h05t10s10d00','barRip1_h10t10s10d00','barRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'terRip1_h05t10s10d00','terRip1_h10t10s10d00','terRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'spreadRip-barRip0'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'barRip0_h10t10s00d00','barRip0_h10t10s02d00','barRip0_h10t10s04d00','barRip0_h10t10s10d00','barRip0_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-barRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'barRip1_h10t10s00d00','barRip1_h10t10s02d00','barRip1_h10t10s04d00','barRip1_h10t10s10d00','barRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-terRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
end



dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat'])
%
%
% marker size, depends on wave-height:
marker_size = 4+4*dat.height;
spread      = dat.spread;

for jj=1:length(spread)

   str = split(dat.run_dirs{jj},'_');
   [~,idx_bathy] = ismember(str{1},possible_bathy);
   idx_spread = find(possible_spreads==spread(jj));
    
   figure(fig0)
   hold on,
   plot( dat.Uscale(jj), dat.Umax(jj),bathy_markers{idx_bathy},'markersize',marker_size(jj),'markerfacecolor',cm(idx_spread,:),'markeredgecolor',cm(idx_spread,:))

   figure(fig1)
   hold on,
   disp('kludging location of Uex')
   idx = find(dat.x>=210,1,'first');
   plot( dat.Uscale(jj), dat.Uex_channel(idx,jj),bathy_markers{idx_bathy},'markersize',marker_size(jj),'markerfacecolor',cm(idx_spread,:),'markeredgecolor',cm(idx_spread,:))

end

end


figure(fig0)
ax = gca;
hold on,
plot([0 1.5], [0 1.5],'--k')
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
axis equal
grid on
xlabel('$\sqrt{-2g\Delta\eta}$ [m/s]','interpreter','latex')
ylabel('$\mathrm{max}(u)$ [m/s]','interpreter','latex')
set(ax,'xlim',[0 1.5], 'ylim',[0 1.5],'ticklabelinterpreter','latex','fontsize',10)
plot(1.18,0.4,'ok','markersize',6), text(1.22,0.4,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
plot(1.18,0.5,'ok','markersiZe',8), text(1.22,0.5,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
plot(1.18,0.6,'ok','markersize',10),text(1.22,0.6,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',8,'location','southeast')

colormap(cm)
cb = colorbar; caxis([0 N])
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
pos = get(ax,'Position');
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'tickdir','out','ticklabelinterpreter','latex','fontsize',10)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

figname = [figDIR,filesep,'Umax_total_vs_Uscale_all.pdf']
exportgraphics(fig0,figname)


figure(fig1)
ax = gca;
hold on,
plot([0 1.5], [0 1.5]/4,'--k')
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
axis equal
grid on
xlabel('$\sqrt{-2g\Delta\eta}$ [m/s]','interpreter','latex')
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(ax,'xlim',[0 1.5], 'ylim',[0 1.5]/4,'ticklabelinterpreter','latex','fontsize',10)
plot(1.18,0.1,'ok','markersize',6), text(1.22,0.1,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
plot(1.18,0.15,'ok','markersiZe',8), text(1.22,0.15,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
plot(1.18,0.2,'ok','markersize',10),text(1.22,0.2,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',8,'location','southeast')

colormap(cm)
cb = colorbar; caxis([0 N])
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
pos = get(ax,'Position');
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'tickdir','out','ticklabelinterpreter','latex','fontsize',10)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

figname = [figDIR,filesep,'Uex_total_vs_Uscale_all.pdf']
exportgraphics(fig1,figname)