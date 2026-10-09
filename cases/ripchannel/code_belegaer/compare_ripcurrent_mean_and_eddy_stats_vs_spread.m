%% there is now a more versitile version of this in ~/git/funwave/cases/ripchannel/code/compare_WaveAvgVelocity_stats_between_runs.m
clear all
close all
addpath('~/git/funwave/code/')
%
%% 1) Estimate the time/space mean cross-shore transport:
%     1.1) time mean transport (in channel vs not)
%     1.2) average over y: ambient region--> 0-1000m and 2000-3000m
%                          ripchannel    --> 1500m +/- L_rc
%% 2) Estimate the time/space eddy cross-shore transport:
%     2.1) Eddy mean transport (in channel vs not)
%     2.2) average over same y ranges
%
%
runBATHY = 'spreadRip'
%
prefixes = {'barRip0','barRip1','terRip1'}; 
suffixes = {'h10t10s00d00','h10t10s02d00','h10t10s04d00','h10t10s10d00','h10t10s20d00'};
mrkrs    = {'o','d','s'};
lbls     = {'Barred $L_y=100$ m','Barred $L_y=50$ m','Terraced $L_y=50$ m'}
%
Ns = length(suffixes);
%
fig0 = figure;% total max & exchange velocities vs velocity scale
fig1 = figure;% eddy  max & exchange velocities vs velocity scale
fig2 = figure;% mean  max & exchange velocities vs velocity scale
fig3 = figure;% mean  exchange in channel and ambient vs spread
fig4 = figure;% eddy  exchange in channel and ambient vs spread
fig5 = figure;% total exchange in channel and ambient vs spread
fig55= figure;% alongshore velocity scale
%
out = struct([]);
for oo = 1:length(prefixes)
runIDs   = cellstr(cell2mat(cat(2,repmat(prefixes(oo),Ns,1),repmat({'_'},Ns,1),suffixes')));
% $$$     runIDs   = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
%
N = length(runIDs);
%
%
dETA          = [];
Tmean_channel = [];
Teddy_channel = [];
T_channel = [];
Umean_channel = [];
Ueddy_channel = [];
U_channel = [];
Tmean_ambient = [];
Teddy_ambient = [];
T_ambient = [];
Umean_ambient = [];
Ueddy_ambient = [];
U_ambient = [];
Umean_rms_channel = [];
Ueddy_rms_channel = [];
Umean_max_channel = [];
Ueddy_max_channel = [];
Utot_max_channel = [];
Vmean         = [];
sig   = {};
%
% $$$ fig5 = figure;
fig6 = figure;
fig7 = figure;
for ii=1:N
    runID = runIDs{ii};
    info = prep_local_ripchannel_info(runBATHY,runID);
    %
    % get value of spread:
    str = split(runID,'_');
    str = split(str{2},{'h','t','s','d'});
    sig{ii} = str{4};
    %
    momFile = dir([info.rootMat,'*',runID,'*MomentumTerms.nc']);
    momFile = [momFile(1).folder,filesep,momFile(1).name]
    if ii==1
        x = ncread(momFile,'x');
        y = ncread(momFile,'y');
        depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
        depFile = [depFile(1).folder,filesep,depFile(1).name];
        dep     = ncread(depFile,'dep');
    end
    %
    % load the "means"
    % Load ETA to get E0 & E1:
    ETA = ncread(momFile,'etamean');
    H   = dep+ETA;
    %
    % Load (U,V) to get T:
    U  = ncread(momFile,'umean');
    V  = ncread(momFile,'vmean');
    %
    % estimate mean transport
    UH = mean(U.*H,3);
    %
    % enforce continuity on mean transport...
    UHtot = mean(UH,1);
    H     = mean(H,3);
    mask  = H>=0.1;
    % convert to an anti-stokes velocity
    Utot  = UHtot./mean(max(H,0.1),1);
    %
    % correct tranport/velocity for stokes transport
    UH = (UH-UHtot).*mask;
    U  = (U-Utot).*mask;
% $$$     Umax = max(U,[],3);
    U  = mean(U,3);
    %
% $$$     % load eddy statistics
% $$$     uFiles = dir([info.rootMat,'*',runID,'*uwavg*.nc']);
% $$$     vFiles = dir([info.rootMat,'*',runID,'*vwavg*.nc']);
% $$$     eFiles = dir([info.rootMat,'*',runID,'*etawavg*.nc']);
% $$$     Nf     = length(uFiles);
% $$$     u = [];
% $$$     v = []; 
% $$$     eta = []; 
% $$$     for jj = 1:Nf
% $$$         u  = cat(3,  u, ncread( [uFiles(jj).folder,filesep,uFiles(jj).name],'uwavg'));
% $$$         eta= cat(3,eta, ncread( [eFiles(jj).folder,filesep,eFiles(jj).name],'etawavg'));
% $$$     end
% $$$     %
% $$$     u = (u-U);
% $$$     disp( 'fixing depth to time-mean... issue in source code' )
% $$$     h  = mean(H,3); %dep+eta;
% $$$     uh = u.*h;
% $$$     %
% $$$     % continuity:
% $$$     uhtot = mean(uh,[1 3],'omitnan');
% $$$     uh    = (uh-uhtot).*mask;
% $$$     utot  = uhtot./mean(max(mean(h, 3),0.1),1);
% $$$     u     = (u-utot).*mask;
% $$$     %
% $$$     % compare continuity with Sxx based estimate:
% $$$     % \bar{u'\eta'}/\bar{h} ~ Hs^2 sqrt(g)/(16 h) ~ 2/3 Sxx/(sqrt(g) h)
% $$$     tmpFig = figure;
% $$$     Sxx = ncread(momFile,'Sxx');
% $$$     uStokes = mean( ( 2./(3*sqrt(9.8)*max(H,0.01)) ).*Sxx, [1 3],'omitnan');
% $$$     plot(x,Utot,'k',x,utot,'--b',x,uStokes,':r')
% $$$     legend({'avg($\langle U \rangle$)','avg($\bar{u}$)','$\approx U_\mathrm{St}$'},'interpreter','latex')
% $$$     xlabel('$x$ [m]','interpreter','latex')
% $$$     ylabel('[m/s]','interpreter','latex')
% $$$     set(gca,'ticklabelinterpreter','latex','tickdir','out','xlim',[50 300])
% $$$     figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'continuity_vs_stokes_',prefixes{oo},'_s',sig{ii},'.pdf'];
% $$$     exportgraphics(tmpFig,figname)
% $$$     close(tmpFig)
    uFiles = dir([info.rootMat,'*',runID,'_velocity_decomposition.nc']);
    Nf     = length(uFiles);    
    u = [];
    for jj = 1:Nf
        u  = cat(3,  u, ncread( [uFiles(jj).folder,filesep,uFiles(jj).name],'Urot'));
    end
    disp( 'fixing depth to time-mean... issue in source code' )
    h  = mean(H,3); %dep+eta;
    uh = u.*h;
    %
    %
    % Estimate alongshore averages:
    iY0 =  y<info.Ly/2 - 5*info.lc | y>info.Ly/2 + 5*info.lc;
    iY1 = (y>info.Ly/2 -   info.lc & y<info.Ly/2 +   info.lc);
    iY2 = (y>info.Ly/2 - 3*info.lc & y<info.Ly/2 + 3*info.lc);
    iX0 = (x>0.5*info.xc & x<info.xc);
    iX1 = find(x>=info.xc,1,'first');
    %
    e0  = mean(ETA(iY0,:,:),[1 3],'omitnan');
    e1  = mean(ETA(iY1,:,:),[1 3],'omitnan');
    %
    Amean_channel =      sum( H(iY1,:,:) , 1, 'omitnan');
    Aeddy_channel =mean( sum( h(iY1,:,:) , 1, 'omitnan') ,3);    
    Amean_ambient =      sum( H(iY0,:,:) , 1, 'omitnan');
    Aeddy_ambient =mean( sum( h(iY0,:,:) , 1, 'omitnan') ,3);    
    %
    Toff       = (UH+uh)>=0 & mask;
    Toff_mean  = UH>=0 & mask;
    Toff_eddy  = uh>=0 & mask;    
    %
    t_channel     = mean( sum((UH(iY1,:,:)+uh(iY1,:,:)).*Toff(iY1,:,:), 1, 'omitnan'), 3);
    tmean_channel =       sum(UH(iY1,:,:).*Toff_mean(iY1,:,:), 1, 'omitnan');
    teddy_channel = mean( sum(uh(iY1,:,:).*Toff_eddy(iY1,:,:), 1, 'omitnan'), 3);
    %
    u_channel     = t_channel./Amean_channel;
    umean_channel = tmean_channel./Amean_channel;
    ueddy_channel = teddy_channel./Aeddy_channel;
    %
    t_ambient     = mean( sum((UH(iY0,:,:)+uh(iY0,:,:)).*Toff(iY0,:,:), 1, 'omitnan'), 3);    
    tmean_ambient =       sum(UH(iY0,:,:).*Toff_mean(iY0,:,:), 1, 'omitnan');
    teddy_ambient = mean( sum(uh(iY0,:,:).*Toff_eddy(iY0,:,:), 1, 'omitnan'), 3);
    %
    u_ambient     = t_ambient./Amean_ambient;
    umean_ambient = tmean_ambient./Amean_ambient;
    ueddy_ambient = teddy_ambient./Aeddy_ambient;
    %
    urms_channel      = rms( U(iY1,:,:)+u(iY1,:,:), [1 3]);
    umean_rms_channel = rms( U(iY1,:,:) , [1 3]);
    ueddy_rms_channel = rms( u(iY1,:,:) , [1 3]);
    urms_ambient      = rms( U(iY0,:,:)+u(iY0,:,:), [1 3]);    
    umean_rms_ambient = rms( U(iY0,:,:) , [1 3]);
    ueddy_rms_ambient = rms( u(iY0,:,:) , [1 3]);
    umax_channel      = max( U(iY1,:,:)+u(iY1,:,:), [], [1 3]);
    umean_max_channel = max( U(iY1,:,:) , [], [1 3]);
    ueddy_max_channel = max( u(iY1,:,:) , [], [1 3]);
    %
    %
    vmean = mean(V(:,iX0,:), [2 3]);
    %
    %% make subplots of each RC transport/max estimate
    iOff = x'>0.7*info.xc;
% $$$     [~,idx_max] = max(umax_channel.*iOff,[],2);
    figure(fig6)
    subplot(N,1,ii)
    imagesc(x,(y(iY2)-info.Ly/2)/info.lc,U(iY2,:)), caxis([-0.5 0.5]), colormap(cmocean('balance'))
% $$$     xline(x(idx_max),'--g')
    xline(x(iX1),'--g')
    yline([-1 1],'--b')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[50 400])
    annotation('textbox','units','normalized','position',[pos(1:2)+[0 0.25].*pos(3:4), 0.5, 0.1],...
               'string',{sprintf('$\\sigma_\\theta=%d^\\circ$',str2num(sig{ii}))},...
               'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',6,'backgroundcolor','none')    
    if ii==1
        title(sprintf('max$(\\langle u\\rangle)$: %s',prefixes{oo}),'interpreter','latex')
        set(gca,'xticklabel',[])
    elseif ii<N
        set(gca,'xticklabel',[])
        if ii==3
            ylabel('$(y-y_0)/L_c$ [m]','interpreter','latex')
        end
    end
    colorbar
    %
    %% make subplots of each RC transport/max estimate
    %    [~,idx_max] = max(ueddy_max_channel.*iOff,[],2);
    figure(fig7)
    subplot(N,1,ii)
    imagesc(x,(y(iY2)-info.Ly/2)/info.lc,max(u(iY2,:,:),[],3)), caxis([-0.5 0.5]), colormap(cmocean('balance'))
% $$$     xline(x(idx_max),'--g')
    xline(x(iX1),'--g')
    yline([-1 1],'--b')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[50 400])
    annotation('textbox','units','normalized','position',[pos(1:2)+[0 0.25].*pos(3:4), 0.5, 0.1],...
               'string',{sprintf('$\\sigma_\\theta=%d^\\circ$',str2num(sig{ii}))},...
               'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',6,'backgroundcolor','none')    
    if ii==1
        title(sprintf('max$(\\bar{u})$: %s',prefixes{oo}),'interpreter','latex')
        set(gca,'xticklabel',[])
    elseif ii<N
        set(gca,'xticklabel',[])
        if ii==3
            ylabel('$(y-y_0)/L_c$ [m]','interpreter','latex')
        end
    end
    colorbar
    %
    %
    dETA         (:,ii) = e1-e0;
    Tmean_channel(:,ii) = tmean_channel;
    Teddy_channel(:,ii) = teddy_channel;
    T_channel    (:,ii) = t_channel;
    Umean_channel(:,ii) = umean_channel;
    Ueddy_channel(:,ii) = ueddy_channel;
    U_channel    (:,ii) = u_channel;
    Urms_channel (:,ii) = urms_channel;    
    Tmean_ambient(:,ii) = tmean_ambient;
    Teddy_ambient(:,ii) = teddy_ambient;
    T_ambient    (:,ii) = t_ambient;
    Umean_ambient(:,ii) = umean_ambient;
    Ueddy_ambient(:,ii) = ueddy_ambient;
    U_ambient    (:,ii) = u_ambient;
    Urms_ambient (:,ii) = urms_ambient;    
    Umean_rms_channel(:,ii) = umean_rms_channel;
    Ueddy_rms_channel(:,ii) = ueddy_rms_channel;    
    Umean_rms_ambient(:,ii) = umean_rms_ambient;
    Ueddy_rms_ambient(:,ii) = ueddy_rms_ambient;    
    Umean_max_channel(:,ii) = umean_max_channel;
    Ueddy_max_channel(:,ii) = ueddy_max_channel;
    Utot_max_channel(:,ii)  = umax_channel;
    Vmean        (:,ii) = vmean;
end
figure(fig6);
xlabel('$x$ [m]','interpreter','latex')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'maximum_channel_mean_speed',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)

figure(fig7);
xlabel('$x$ [m]','interpreter','latex')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'maximum_channel_eddy_speed',prefixes{oo},'.pdf'];
exportgraphics(fig7,figname)

close(fig6)
close(fig7)

cm = cmocean('thermal',N+1);
cm = cm(1:N,:);


fig8 = figure;
colororder(cm);
ax1 = subplot(2,1,1);
p1 = plot(x,U_channel,'-');
xline(50,'--r')
xline(x(iX1),'--g')
title(sprintf('Channel: %s',prefixes{oo}))
ylabel('$U_\mathrm{tot,ex}$ [m/s]~~~~~~~~~~~~~~','interpreter','latex')
legend(p1([1]),{'channel'},'interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.2])
ax2 = subplot(2,1,2);
p2 = plot(x,U_ambient,'-');
xline(50,'--r')
xline(x(iX1),'--g')
legend(p2([1]),{'ambient'},'interpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.2])
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax1.Position(3) = ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_total_channel_vs_ambient_',prefixes{oo},'.pdf'];
exportgraphics(fig8,figname)
close(fig8)


fig9 = figure;
colororder(cm);
ax1 = subplot(2,1,1);
p1 = plot(x,Umean_channel,'-');
xline(50,'--r')
xline(x(iX1),'--g')
title(sprintf('Channel: %s',prefixes{oo}))
ylabel('$U_\mathrm{ex}$ [m/s]~~~~~~~~~~~~~~~~~','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.2])
legend(p1([1]),{'$\langle{u}\rangle$'},'interpreter','latex')
ax2 = subplot(2,1,2);
p2 = plot(x,Ueddy_channel,'-');
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.2])
legend(p2([1]),{'$\bar{u}$'},'interpreter','latex')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax1.Position(3) = ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_channel_mean_vs_eddy_',prefixes{oo},'.pdf'];
exportgraphics(fig9,figname)
close(fig9)
%
fig10 = figure;
colororder(cm);
ax1 = subplot(2,1,1);
p1 = plot(x,Umean_ambient,'-');
xline(50,'--r')
xline(x(iX1),'--g')
ylabel('$U_\mathrm{ex}$ [m/s]~~~~~~~~~~~~~~~~~','interpreter','latex')
title(sprintf('Ambient: %s',prefixes{oo}))
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.1])
legend(p1([1]),{'$\langle u\rangle$'},'interpreter','latex')
ax2 = subplot(2,1,2);
p2 = plot(x,Ueddy_ambient,'-');
xline(50,'--r')
xline(x(iX1),'--g')
xlabel('$x$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.1])
legend(p2([1]),{'$\bar{u}$'},'interpreter','latex')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax1.Position(3) = ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_ambient_mean_vs_eddy_',prefixes{oo},'.pdf'];
exportgraphics(fig10,figname)
close(fig10)
%
fig11 = figure;
colororder(cm);
% $$$ iOff = x>0.7*info.xc;
% $$$ [~,idx_max1] = max(Utot_max_channel.*iOff,[],1);
% $$$ idx_max2 = sub2ind(size(Ueddy_max_channel),idx_max1,1:N);
ax0 = subplot(3,1,1);
p0 = plot(x,Utot_max_channel,'-'); hold on,
xline(50,'--r')
xline(x(iX1),'--g')
ax1 = subplot(3,1,2);
p1 = plot(x,Umean_max_channel,'-'); hold on,
xline(50,'--r')
xline(x(iX1),'--g')
ax2 = subplot(3,1,3);
p2 = plot(x,Ueddy_max_channel,'-'); hold on,
xline(50,'--r')
xline(x(iX1),'--g')
%
% $$$ for kk=1:N
% $$$     axes(ax0)
% $$$     plot(x(iX1),Umax_channel(iX1,kk),'x','color',cm(kk,:))
% $$$     axes(ax1)
% $$$     plot(x(iX1),Umean_max_channel(iX1,kk),'x','color',cm(kk,:))
% $$$     axes(ax2)
% $$$     plot(x(iX1),Ueddy_max_channel(iX1,kk),'x','color',cm(kk,:))        
% $$$ % $$$     plot(x(idx_max1(kk)),Umean_max_channel(idx_max2(kk)),'x',x(idx_max1(kk)),Ueddy_max_channel(idx_max2(kk)),'x','color',cm(kk,:))
% $$$ end
%
axes(ax0)
set(gca,'ticklabelinterpreter','latex','ylim',[0 max(1,max(Ueddy_max_channel(:)))])
title(sprintf('%s',prefixes{oo}))
legend(p0([1]),{'$\langle u\rangle+\bar{u}$'},'interpreter','latex')
axes(ax1)
ylabel('$U_\mathrm{max}$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 max(1,max(Ueddy_max_channel(:)))])
title(sprintf('%s',prefixes{oo}))
legend(p1([1]),{'$\langle u\rangle$'},'interpreter','latex')
axes(ax2)
xlabel('$x$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 max(1,max(Ueddy_max_channel(:)))])
legend(p2([1]),{'$\bar{u}$'},'interpreter','latex')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax0.Position(3)=ax2.Position(3);
ax1.Position(3)=ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Umax_channel_mean_vs_eddy_',prefixes{oo},'.pdf'];
exportgraphics(fig11,figname)
close(fig11)
%
%
fig12 = figure;
colororder(cm);
ax1 = subplot(2,1,1);
p1 = plot(x,Umean_rms_channel,'-');
xline(50,'--r')
xline(x(iX1),'--g')
ylabel('$U_\mathrm{rms}$ [m/s]~~~~~~~~~~~~~~~~~','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.5])
legend(p1([1]),{'$\langle u\rangle$'},'interpreter','latex')
title(sprintf('Channel: %s',prefixes{oo}))
ax2 = subplot(2,1,2);
p2 = plot(x,Ueddy_rms_channel,'-');
xline(50,'--r')
xline(x(iX1),'--g')
xlabel('$x$ [m]','interpreter','latex')
legend(p2([1]),{'$\bar{u}$'},'interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.5])
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax1.Position(3)=ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'rms_ripchannel_speed_mean_vs_eddy_',prefixes{oo},'.pdf'];
exportgraphics(fig12,figname)
close(fig12)

fig13 = figure;
colororder(cm);
ax1 = subplot(2,1,1);
p1 = plot(x,Umean_rms_ambient,'-');
xline(50,'--r')
xline(x(iX1),'--g')
ylabel('$U_\mathrm{rms}$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.5])
title(sprintf('Ambient: %s',prefixes{oo}))
legend(p1([1]),{'$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex')
ax2 = subplot(2,1,2);
p2  = plot(x,Ueddy_rms_ambient,'-');
xline(50,'--r')
xline(x(iX1),'--g')
xlabel('$x$ [m]','interpreter','latex')
legend(p2([1]),{'$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[0 0.5])
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
ax1.Position(3)=ax2.Position(3);
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'rms_ambient_speed_mean_vs_eddy_',prefixes{oo},'.pdf'];
exportgraphics(fig13,figname)
close(fig13)

fig14 = figure;
colororder(cm);
p1 = plot((y-info.Ly/2)/info.lc,Vmean,'-');
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
ylabel('$V$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',[-0.25 0.25],'xlim',[-1 1]*7)
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
title(sprintf('Mid-SZ: %s',prefixes{oo}))
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'mean_alongshore_current_speed_',prefixes{oo},'.pdf'];
exportgraphics(fig14,figname)
close(fig14)


    out(oo).dETA         = dETA;
    out(oo).Tmean_channel=Tmean_channel;
    out(oo).Teddy_channel=Teddy_channel;
    out(oo).T_channel    =T_channel;    
    out(oo).Umean_channel=Umean_channel;
    out(oo).Ueddy_channel=Ueddy_channel;
    out(oo).U_channel=U_channel;    
    out(oo).Tmean_ambient=Tmean_ambient;
    out(oo).Teddy_ambient=Teddy_ambient;
    out(oo).T_ambient=T_ambient;    
    out(oo).Umean_ambient=Umean_ambient;
    out(oo).Ueddy_ambient=Ueddy_ambient;
    out(oo).U_ambient=U_ambient;    
    out(oo).Urms_channel=Urms_channel;
    out(oo).Umean_rms_channel=Umean_rms_channel;
    out(oo).Ueddy_rms_channel=Ueddy_rms_channel;
    out(oo).Urms_ambient=Urms_ambient;              
    out(oo).Umean_rms_ambient=Umean_rms_ambient;
    out(oo).Ueddy_rms_ambient=Ueddy_rms_ambient;
    out(oo).Umean_max_channel=Umean_max_channel;
    out(oo).Ueddy_max_channel=Ueddy_max_channel;

    % maximum ambient exchange velocity
    iOff = (x>=(info.xc-info.dx) & x<=(info.xc+info.dx));
    Uex_ambient      = max( U_ambient.*iOff, [], 1);    
    Uex_mean_ambient = max( Umean_ambient.*iOff, [], 1);
    Uex_eddy_ambient = max( Ueddy_ambient.*iOff, [], 1);    
    
    % rip-current transport velocity scale
    Uex_channel      = max( U_channel.*iOff, [], 1);
    Uex_mean_channel = max( Umean_channel.*iOff, [], 1);
    Uex_eddy_channel = max( Ueddy_channel.*iOff, [], 1);

    % get fastest offshore speed
    [Umax, idx_max] = max(Utot_max_channel.*iOff,[],1);
    idx_max2       = sub2ind(size(Ueddy_max_channel),idx_max,1:N);
    Umax_mean      = Umean_max_channel(idx_max2);
    Umax_eddy      = Ueddy_max_channel(idx_max2);
    
    % U-scale:
    g = 9.8;
    Uscale = sqrt(2*g*abs(mean(dETA(iX0,:))));
    %
    %
    % V-scale:
    iYp = y>info.Ly/2 & y<info.Ly/2+5*info.lc;
    iYm = y<info.Ly/2 & y>info.Ly/2-5*info.lc;
    Vscale = 0.5*(-min(Vmean(iYp,:)) + max(Vmean(iYm,:)));

    % log values
    Uex_mean_ambient_log(:,oo) = Uex_mean_ambient;
    Uex_eddy_ambient_log(:,oo) = Uex_eddy_ambient;
    Uex_ambient_log     (:,oo) = Uex_ambient;    
    Uex_mean_channel_log(:,oo) = Uex_mean_channel;
    Uex_eddy_channel_log(:,oo) = Uex_eddy_channel;
    Uex_channel_log     (:,oo) = Uex_channel;    
    Xmax_log     (:,oo)        = x(idx_max);
    Umax_mean_log(:,oo)        = Umax_mean;
    Umax_eddy_log(:,oo)        = Umax_eddy;
    Umax_log     (:,oo)        = Umax;
    Uscale_log   (:,oo)        =Uscale;
    Vscale_log   (:,oo)        =Vscale;
    
    %% plot total maximum and exchange ripchannel speed vs. velocity scale
    figure(fig0),
    if oo==1
        plot([0 1],[0 1],'--k')
        colororder(cm)
    end
    hold on,
    for ii=1:N
        plot(Uscale(ii),Umax(ii),mrkrs{oo},'color',cm(ii,:),'markerfacecolor',cm(ii,:))
    end
    
    for ii=1:N
        plot(Uscale(ii),Uex_channel(ii),mrkrs{oo},'color',cm(ii,:))
    end


    %% plot eddy maximum and exchange ripchannel speed vs. velocity scale
    figure(fig1),
    if oo==1
        plot([0 1],[0 1],'--k')
        colororder(cm)
    end
    hold on,
    for ii=1:N
        plot(Uscale(ii),Umax_eddy(ii),mrkrs{oo},'color',cm(ii,:),'markerfacecolor',cm(ii,:))
    end
    
    for ii=1:N
        plot(Uscale(ii),Uex_eddy_channel(ii),mrkrs{oo},'color',cm(ii,:))
    end


    %% plot mean maximum and exchange ripchannel speed vs. velocity scale
    figure(fig2),
    if oo==1
        plot([0 1],[0 1],'--k')
        colororder(cm)
    end
    hold on,
    for ii=1:N
        plot(Uscale(ii),Umax_mean(ii),mrkrs{oo},'color',cm(ii,:),'markerfacecolor',cm(ii,:))
    end
    
    for ii=1:N
        plot(Uscale(ii),Uex_mean_channel(ii),mrkrs{oo},'color',cm(ii,:))
    end
    
    %% plot mean Uex in channel vs ambient
    figure(fig3), hold on,
    if oo==1
        colororder(cm)
    end
    spread = str2num(cell2mat(sig'));    
    plot(spread,Uex_mean_channel,mrkrs{oo},'color','k','markerfacecolor','k')
    plot(spread,Uex_mean_ambient,mrkrs{oo},'color','r','markerfacecolor','r')
    

    %% plot eddy Uex in channel vs ambient
    figure(fig4), hold on,
    if oo==1
        colororder(cm)
    end
    plot(spread,Uex_eddy_channel,mrkrs{oo},'color','k','markerfacecolor','k')
    plot(spread,Uex_eddy_ambient,mrkrs{oo},'color','r','markerfacecolor','r')

    figure(fig5), hold on,
    if oo==1
        colororder(cm)
    end
    plot(spread,Uex_channel,mrkrs{oo},'color','k','markerfacecolor','k')
    plot(spread,Uex_ambient,mrkrs{oo},'color','r','markerfacecolor','r')

    %% plot eddy Uex in channel vs ambient
    figure(fig55), hold on,
    if oo==1
        colororder(cm)
    end
    plot(spread,Vscale,mrkrs{oo},'color','k','markerfacecolor','k')
    plot(spread,Uscale,mrkrs{oo},'color','r','markerfacecolor','r')

end

figure(fig0)
xlabel('$\sqrt{2g\Delta\eta_y}$ [m/s]','interpreter','latex')
ylabel('$(\langle u\rangle + \bar{u})$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
title('Maximum (filled) and Exchange (open) Velocity Scales','interpreter','latex','fontsize',10)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
legend(gca().Children([5*N, 3*N, N]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Umax_and_Uex_total_ripcurrent_speed_vs_velocity_scale.pdf'];
exportgraphics(fig0,figname)

figure(fig1)
xlabel('$\sqrt{2g\Delta\eta_y}$ [m/s]','interpreter','latex')
ylabel('$\bar{u}$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
title('Maximum (filled) and Exchange (open) Velocity Scales','interpreter','latex','fontsize',10)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
legend(gca().Children([5*N, 3*N, N]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Umax_and_Uex_eddy_ripcurrent_speed_vs_velocity_scale.pdf'];
exportgraphics(fig1,figname)

figure(fig2)
xlabel('$\sqrt{2g\Delta\eta_y}$ [m/s]','interpreter','latex')
ylabel('$\langle u\rangle$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
title('Maximum (filled) and Exchange (open) Velocity Scales','interpreter','latex','fontsize',10)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
legend(gca().Children([5*N, 3*N, N]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Umax_and_Uex_mean_ripcurrent_speed_vs_velocity_scale.pdf'];
exportgraphics(fig2,figname)


figure(fig3)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$\langle{u}\rangle$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 0.25])
title('Channel (black) and Ambient (red) Exchange Velocity Scales','interpreter','latex','fontsize',10)
legend(gca().Children([5, 3, 1]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_mean_channel_and_ambient_vs_spread.pdf'];
exportgraphics(fig3,figname)

figure(fig4)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$\bar{u}$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 0.25])
title('Channel (black) and Ambient (red) Exchange Velocity Scales','interpreter','latex','fontsize',10)
legend(gca().Children([5, 3, 1]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_eddy_channel_and_ambient_vs_spread.pdf'];
exportgraphics(fig4,figname)

figure(fig5)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$\langle{u}\rangle + \bar{u}$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 0.25])
title('Channel (black) and Ambient (red) Exchange Velocity Scales','interpreter','latex','fontsize',10)
legend(gca().Children([5, 3, 1]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Uex_total_channel_and_ambient_vs_spread.pdf'];
exportgraphics(fig5,figname)

figure(fig55)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('[m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 0.5])
title('V (black) and $\sqrt{2g\Delta\langle{\eta}\rangle_y}$ (red) Velocity Scales','interpreter','latex','fontsize',10)
legend(gca().Children([5, 3, 1]),lbls,'interpreter','latex','autoupdate','off','fontsize',10)

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'Vscale_and_Uscale_vs_spread.pdf'];
exportgraphics(fig55,figname)


% close all
save(['/data2/ripchannel/',runBATHY,'/mat_data/mean_vs_eddy_rip_transport_and_speed.mat'],'out','Uex_mean_ambient_log','Uex_eddy_ambient_log','Uex_ambient_log','Uex_mean_channel_log','Uex_eddy_channel_log','Uex_channel_log','Xmax_log','Umax_mean_log','Umax_eddy_log','Umax_log','Uscale_log','Vscale_log')