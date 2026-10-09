%% create compendium plot of exchange and enstrophy for a set of runs:
clear all
close all
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'

%% size for height
possible_heights = [0.5 1.0 1.5];
marker_sizes     = [2 4 6];
cm_height = cmocean('thermal',4);
cm_height = cm_height(1:end-1,:);

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

%% 2) velocity scales based on high-rip evolution:
BATHY = {};
C_MOD = {};% scales based on modeled
C_EST = {};% scales based on offshore
C_EST_EX = {};% scales based on zero-spread for exchange
C_EST_RMS = {};% scales based on zero-spread for exchange


% figure parameters
xm = 2.5;
ym = 2.5;
pw = 4;
ph = 4;
ag = 0.5;
ppos1  = [xm       ym         pw ph];
cbpos  = [xm+pw+ag ym       ag ph/2];
ps     = [2*xm+pw+6*ag  2*ym+(ag+ph)];

fig0 = figure('units','centimeters');% < max(u)  > vs Uscale
fig0.Position(3:4)=ps;
set(fig0,'papersize',ps,'paperposition',[0 0 ps])
ax0 = axes('units','centimeters','position',ppos1);
%
fig1 = figure('units','centimeters');% < max(Uex)> vs Uscale
fig1.Position(3:4)=ps;
set(fig1,'papersize',ps,'paperposition',[0 0 ps])
ax1 = axes('units','centimeters','position',ppos1);

fig2 = figure;% < max(u)  >/Uscale_offshore vs Spread/Hs
fig3 = figure;% < max(u)  >/Uscale_observed vs Spread/Hs
fig4 = figure;% Uscale_observed vs Uscale_offshore
fig5 = figure;% Uex vs Uscale_offshore
fig6 = figure;% Urms vs Uscale_offshore
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



dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat']);
%
%
% marker size, depends on wave-height:
height      = dat.height;
marker_size = 1+4*height;
spread      = dat.spread;

split_name = split(NAME,'-')
if ismember('highRip',split_name) & ismember('s00',split_name)
    %% for safety... don't divide by zero:
    tmp = dat.Umax'\dat.Uscale';
    BATHY = cat(1,BATHY,split_name{2});
    C_MOD = cat(1,C_MOD,1/tmp);
    clear Uscale0
    for jj=1:length(dat.height)
        %%
        %% Plot Moulton's scaling:
        g = 9.8;
        gam = 0.35;
        %% Kludge to find bar-crest:
        icrest = find(dat.x>=dat.bar_location(jj),1,'first');
        hcrest = 0.03*(dat.x(icrest)-50)-dat.bar_amplitude(jj).*exp(-0.5.*((dat.x(icrest)-dat.bar_location(jj))./dat.bar_width(jj)).^2);
        %
        %% attempt 1 to estimate Hbr:
        k0 = wavenumber_FunwaveTVD(2*pi/dat.period(jj),9);
        disper0 = sqrt(g.*k0.*tanh(k0.*9));
        dwdk_full0 = 0.5*(g*tanh(k0.*9) + g*k0.*9.*sech(k0.*9).^2)./disper0;
        Hbr0 = nthroot( dat.height(jj).^4.*dwdk_full0.^2.*(gam/g), 5);
% $$$         %% attempt 2 to estimate Hbr
% $$$         k  = wavenumber_FunwaveTVD(2*pi/dat.period(jj),dat.h0);
% $$$         disper = sqrt(g.*k.*tanh(k.*dat.h0));
% $$$         dwdk_full = 0.5*(g*tanh(k.*dat.h0) + g*k.*dat.h0.*sech(k.*dat.h0).^2)./disper;
% $$$         Hshoal = real(sqrt( dat.height(jj).^2.*dwdk_full(end)./dwdk_full ));
% $$$         Hmax   = 0.88./k.*tanh(gam*k.*dat.h0/0.88);
% $$$         ibrk   = find(Hshoal>=Hmax,1,'last');
% $$$         Hbr = Hmax(ibrk);
        dETA_moulton0 = -gam.^2/16.*(cosd(dat.direction(jj)).^2 + 0.5).*max(Hbr0/gam-hcrest,0);
% $$$         dETA_moulton1 = -gam.^2/16.*(cosd(dat.direction(jj)).^2 + 0.5).*max(Hbr/gam-hcrest ,0);       
        Uscale0(jj) = sqrt(-2*g*dETA_moulton0);
% $$$         Uscale1(jj) = sqrt(-2*g*dETA_moulton1);       
%
%
%
        %% compile exchange velocities to scale zero-spread
        iX  = find(dat.x>=dat.Xke_max(jj),1,'first');
        Uex(jj) = dat.Uex_channel(iX,jj);
        
        %%
    end
    tmp   = Uscale0/dat.Umax;
    C_EST = cat(1,C_EST,1/tmp);
    
    tmp = Uscale0/dat.Uke_channel_max;
    C_EST_RMS = cat(1,C_EST_RMS,1/tmp);

    tmp   = Uscale0/Uex;
    C_EST_EX  = cat(1,C_EST_EX,1/tmp);
end


for jj=1:length(spread)

   str = split(dat.run_dirs{jj},'_');
   [~,idx_bathy] = ismember(str{1},possible_bathy);
   idx_spread = find(possible_spreads==spread(jj));
   idx_height = find(possible_heights==height(jj));   

   [~,idx_Uscale] = ismember(str{1},BATHY);

   %% Plot Moulton's scaling:
   g = 9.8;
   gam = 0.35;
   %% Kludge to find bar-crest:
   icrest = find(dat.x>=dat.bar_location(jj),1,'first');
   hcrest = 0.03*(dat.x(icrest)-50)-dat.bar_amplitude(jj).*exp(-0.5.*((dat.x(icrest)-dat.bar_location(jj))./dat.bar_width(jj)).^2);
   %% attempt 1 to estimate Hbr:
   k0 = wavenumber_FunwaveTVD(2*pi/dat.period(jj),9);
   disper0 = sqrt(g.*k0.*tanh(k0.*9));
   dwdk_full0 = 0.5*(g*tanh(k0.*9) + g*k0.*9.*sech(k0.*9).^2)./disper0;
   Hbr0 = nthroot( dat.height(jj).^4.*dwdk_full0.^2.*(gam/g), 5);
   dETA_moulton0 = -gam.^2/16.*(cosd(dat.direction(jj)).^2 + 0.5).*max(Hbr0/gam-hcrest,0);
   Uscale0 = sqrt(-2*g*dETA_moulton0);
   
   figure(fig0)
   hold on,
   plot( Uscale0.*C_EST{idx_Uscale}, dat.Umax(jj),bathy_markers{idx_bathy},'markersize',marker_size(jj),'markeredgecolor',cm(idx_spread,:))%'markerfacecolor',cm(idx_spread,:),

   figure(fig1)
   hold on,
   disp('kludging location of Uex')
   % UexMax = max(dat.Uex_channel(:,jj).*(dat.x>110));
   %   iX = find(x>=dat.bar_location(jj),1,'first');
   iX  = find(dat.x>=dat.Xke_max(jj),1,'first');
   Uex = dat.Uex_channel(iX,jj);
   plot( Uscale0.*C_EST_EX{idx_Uscale}, Uex,bathy_markers{idx_bathy},'markersize',marker_size(jj),'markeredgecolor',cm(idx_spread,:))% ,'markerfacecolor',cm(idx_spread,:)


   
   figure(fig2)
   hold on,
   plot(dat.spread(jj),dat.Umax(jj)./(C_MOD{idx_Uscale}*dat.Uscale(jj)),bathy_markers{idx_bathy},'markerfacecolor',cm_height(idx_height,:),'markeredgecolor',cm_height(idx_height,:))
   yline(1,'--r')
   
   figure(fig3)
   hold on,
   plot(dat.spread(jj),dat.Umax(jj)./(C_EST{idx_Uscale}*Uscale0),bathy_markers{idx_bathy},'markerfacecolor',cm_height(idx_height,:),'markeredgecolor',cm_height(idx_height,:))
   yline(1,'--r')

   figure(fig4)
   hold on,
   plot((C_EST{idx_Uscale}*Uscale0),(C_MOD{idx_Uscale}*dat.Uscale(jj)),bathy_markers{idx_bathy},'markersize',marker_size(jj),'markerfacecolor',cm(idx_spread,:),'markeredgecolor',cm(idx_spread,:))

   
   figure(fig5)
   hold on,
   plot(dat.spread(jj),Uex./(C_EST_EX{idx_Uscale}*Uscale0),bathy_markers{idx_bathy},'markerfacecolor',cm_height(idx_height,:),'markeredgecolor',cm_height(idx_height,:))


   figure(fig6)
   hold on,
   plot(dat.spread(jj),dat.Uke_channel_max(jj)./(C_EST_RMS{idx_Uscale}*Uscale0),bathy_markers{idx_bathy},'markerfacecolor',cm_height(idx_height,:),'markeredgecolor',cm_height(idx_height,:))

end


% $$$ fig = figure;
% $$$ plot(spread,dat.Uscale,'xk',spread,Uscale0,'sb',spread,Uscale1,'or')

end


figure(fig0)
hold on,
plot([0 1.5], [0 1.5],'--k')
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
axis equal
grid on
xlabel(ax0,'$\mathcal{U}$ [m/s]','interpreter','latex')
ylabel(ax0,'$u_\mathrm{rip}$ [m/s]','interpreter','latex')
set(ax0,'xlim',[0 1], 'ylim',[0 1],'ticklabelinterpreter','latex','fontsize',10)
plot(0.6,0.3,'ok','markersize'  ,marker_sizes(1)), text(0.65,0.3,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',5)
plot(0.6,0.375,'ok','markersiZe',marker_sizes(2)), text(0.65,0.375,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',5)
plot(0.6,0.45,'ok','markersize' ,marker_sizes(3)),text(0.65,0.45,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',5)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',5,'location','southeast')
% $$$ 
% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')
% $$$ pos = get(ax,'Position');
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
% $$$ cb.Position(4) = 0.5*cb.Position(4);
% $$$ ax.Position = pos;
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Umax_total_vs_Uscale_all.pdf']
exportgraphics(fig0,figname)


figure(fig1)
hold on,
plot([0 1.5], [0 1.5],'--k')
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
grid on
xlabel('$c_\mathrm{ex}\mathcal{U}$ [m/s]','interpreter','latex')
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(ax1,'xlim',[0 0.2], 'ylim',[0 0.2],'ticklabelinterpreter','latex','fontsize',10,'dataaspectratio',[1 1 1],'plotboxaspectratio',[1 1 1])
% $$$ plot(0.125,0.06  ,'ok','markersize',6) ,  text(0.13,0.06, sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
% $$$ plot(0.125,0.08  ,'ok','markersiZe',8) ,  text(0.13,0.08,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
% $$$ plot(0.125,0.1   ,'ok','markersize',10),  text(0.13,0.1 ,  sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

% $$$ legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',8,'location','southeast')

% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')
% $$$ pos = get(ax,'Position');
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
% $$$ cb.Position(4) = 0.5*cb.Position(4);
% $$$ ax.Position = pos;

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Uex_total_vs_Uscale_all.pdf']
exportgraphics(fig1,figname)



figure(fig2)
ax = gca;
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
grid on
xlabel('$\sigma_\theta$','interpreter','latex')
ylabel('$u_\mathrm{rip}/(c_\mathrm{rip}\mathcal{U}_\mathrm{mod})$ []','interpreter','latex')
set(ax,'xlim',[0 20], 'ylim',[0 1.4],'ticklabelinterpreter','latex','fontsize',15)
% $$$ plot(15,0.3,'ok','markersize',6), text(16,0.3,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.4,'ok','markersiZe',8), text(16,0.4,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.5,'ok','markersize',10),text(16,0.5,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',10,'location','southeast')

colormap(cm_height)
M = length(possible_heights);
cb = colorbar; caxis([0 M])
ylabel(cb,'$H_s$ [m]','interpreter','latex')
pos = get(ax,'Position');
set(cb,'ylim',[0 M],'ytick',(1:M)-0.5,'yticklabel',num2str(possible_heights'),'tickdir','out','ticklabelinterpreter','latex','fontsize',10)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

figname = [figDIR,filesep,'Umax_total_divided_by_Uscale_MOD.pdf']
exportgraphics(fig2,figname)


figure(fig3)
ax = gca;
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
grid on
xlabel('$\sigma_\theta$','interpreter','latex')
ylabel('$u_\mathrm{rip}/\mathcal{U}$ []','interpreter','latex')
set(ax,'xlim',[0 20], 'ylim',[0 1.4],'ticklabelinterpreter','latex','fontsize',15)
% $$$ plot(15,0.3,'ok','markersize',6), text(16,0.3,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.4,'ok','markersiZe',8), text(16,0.4,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.5,'ok','markersize',10),text(16,0.5,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',10,'location','southeast')

colormap(cm_height)
cb = colorbar; caxis([0 M])
pos = get(ax,'Position');
ylabel(cb,'$H_s$ [m]','interpreter','latex')
set(cb,'ylim',[0 M],'ytick',(1:M)-0.5,'yticklabel',num2str(possible_heights'),'tickdir','out','ticklabelinterpreter','latex','fontsize',10)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

% figname = [figDIR,filesep,'Umax_total_divided_by_Uscale_EST_no_coeff.pdf']
figname = [figDIR,filesep,'Umax_total_divided_by_Uscale_EST.pdf']
exportgraphics(fig3,figname)



figure(fig4)
ax = gca;
hold on,
plot([0 1.5], [0 1.5],'--k')
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
axis equal
grid on
xlabel('$\mathcal{U}$ [m/s]','interpreter','latex')
ylabel('$\mathcal{U}_\mathrm{mod}$ [m/s]','interpreter','latex')
set(ax,'xlim',[0 1.5], 'ylim',[0 1.5],'ticklabelinterpreter','latex','fontsize',15)
plot(1.18,0.4,'ok','markersize',6), text(1.22,0.4,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
plot(1.18,0.5,'ok','markersiZe',8), text(1.22,0.5,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
plot(1.18,0.6,'ok','markersize',10),text(1.22,0.6,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',10,'location','southeast')

colormap(cm)
cb = colorbar; caxis([0 N])
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
pos = get(ax,'Position');
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(possible_spreads'),'tickdir','out','ticklabelinterpreter','latex','fontsize',12)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

figname = [figDIR,filesep,'Umod_vs_Uest_all.pdf']
exportgraphics(fig4,figname)





figure(fig5)
ax = gca;
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
grid on
xlabel('$\sigma_\theta$','interpreter','latex')
ylabel('$U_\mathrm{ex}/(c_\mathrm{ex}\mathcal{U})$ []','interpreter','latex')
%ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(ax,'xlim',[0 20], 'ylim',[0 2.5],'ticklabelinterpreter','latex','fontsize',15)
% $$$ plot(15,0.3,'ok','markersize',6), text(16,0.3,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.4,'ok','markersiZe',8), text(16,0.4,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.5,'ok','markersize',10),text(16,0.5,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',10,'location','southeast')

colormap(cm_height)
cb = colorbar; caxis([0 M])
pos = get(ax,'Position');
ylabel(cb,'$H_s$ [m]','interpreter','latex')
set(cb,'ylim',[0 M],'ytick',(1:M)-0.5,'yticklabel',num2str(possible_heights'),'tickdir','out','ticklabelinterpreter','latex','fontsize',12)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

figname = [figDIR,filesep,'Uex_total_divided_by_Uscale_vs_spread.pdf']
% figname = [figDIR,filesep,'Uex_total_vs_spread.pdf']
exportgraphics(fig5,figname)



figure(fig6)
ax = gca;
p0 = plot(-999,-999,'ok',-999,-999,'dk',-999,-999,'sk');
grid on
xlabel('$\sigma_\theta$','interpreter','latex')
ylabel('$\mathrm{rms}(u)/c_\mathrm{rms}\mathcal{U}$ []','interpreter','latex')
%ylabel('rms$(u)$ [m/s]','interpreter','latex')
set(ax,'xlim',[0 20], 'ylim',[0 1.5],'ticklabelinterpreter','latex','fontsize',15)
% $$$ plot(15,0.25,'ok','markersize',6), text(16,0.25,sprintf('$H_s=%1.1f$ m',possible_heights(1)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.3,'ok','markersiZe',8), text(16,0.3,sprintf('$H_s=%1.1f$ m',possible_heights(2)),'interpreter','latex','fontsize',8)
% $$$ plot(15,0.35,'ok','markersize',10),text(16,0.35,sprintf('$H_s=%1.1f$ m',possible_heights(3)),'interpreter','latex','fontsize',8)

legend(p0,{'Barred 100-m','Barred 50-m','Terraced 50-m'},'interpreter','latex','fontsize',10,'location','southeast')

colormap(cm_height)
cb = colorbar; caxis([0 M])
pos = get(ax,'Position');
ylabel(cb,'$H_s$ [m]','interpreter','latex')
set(cb,'ylim',[0 M],'ytick',(1:M)-0.5,'yticklabel',num2str(possible_heights'),'tickdir','out','ticklabelinterpreter','latex','fontsize',12)
cb.Position(4) = 0.5*cb.Position(4);
ax.Position = pos;

% figname = [figDIR,filesep,'Urms_total_vs_spread.pdf']
figname = [figDIR,filesep,'Urms_total_divided_by_Uscale_vs_spread.pdf']
exportgraphics(fig6,figname)

disp('c_EST:')
disp(C_EST)
disp('c_MOD:')
disp(C_MOD)
disp('C_EST_RMS')
C_EST_RMS
disp('C_EST_EX')
C_EST_EX