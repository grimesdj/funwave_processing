%% create compendium plot of exchange and enstrophy for a set of runs:
clear all
close all
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'

%% 1) need a list of run-directories
NAMES    = {'uniRip-ter2D','uniRip-bar2D','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10','spreadRip-barRip0','spreadRip-barRip1','spreadRip-terRip1'};
for nn = 1:length(NAMES)
    clearvars -except rootDIR figDIR outDIR NAMES nn
    
NAME     = NAMES{nn};
switch NAME
  case 'uniRip-ter2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'}
    run_dirs = {'ter2D_h10t10s02d00','ter2D_h10t10s04d00','ter2D_h10t10s10d00','ter2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'uniRip-bar2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'}    
    run_dirs = {'bar2D_h10t10s02d00','bar2D_h10t10s04d00','bar2D_h10t10s10d00','bar2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'highRip-barRip0-s00'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'barRip0_h05t10s00d00','barRip0_h10t10s00d00','barRip0_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'barRip1_h05t10s00d00','barRip1_h10t10s00d00','barRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'terRip1_h05t10s00d00','terRip1_h10t10s00d00','terRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip0-s10'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'barRip0_h05t10s10d00','barRip0_h10t10s10d00','barRip0_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'barRip1_h05t10s10d00','barRip1_h10t10s10d00','barRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'}    
    run_dirs = {'terRip1_h05t10s10d00','terRip1_h10t10s10d00','terRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'spreadRip-barRip0'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'}
    run_dirs = {'barRip0_h10t10s00d00','barRip0_h10t10s02d00','barRip0_h10t10s04d00','barRip0_h10t10s10d00','barRip0_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-barRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'}
    run_dirs = {'barRip1_h10t10s00d00','barRip1_h10t10s02d00','barRip1_h10t10s04d00','barRip1_h10t10s10d00','barRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-terRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'}
    run_dirs = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
end


%% 2) loop over the list, loading and allocating Uex, ENS, etc.
% wave and grid parameters
height = [];
spread = [];
period = [];
direction = [];
bar_width   = [];
bar_amplitude    = [];
bar_location = [];
channel_length   = [];
channel_amplitude_ratio  = [];
%
N = length(runIDs);
fig0  = figure;
fig00 = figure;
for ii = 1:N
% ii=1
runID  = runIDs{ii};
runDIR = run_dirs{ii};
info   = prep_belegaer_ripchannel_info(runID,runDIR);
%
waves = split(info.runName,'_');
waves  = split(waves{2},{'h','t','s','d'});
%
height = cat(2,height,str2num(waves{2})/10);
period = cat(2,period,str2num(waves{3}));
spread = cat(2,spread,str2num(waves{4}));
direction = cat(2,direction,str2num(waves{5}));
bar_width   = cat(2,bar_width,info.wc);
bar_amplitude    = cat(2,bar_amplitude,info.ac);
bar_location = cat(2,bar_location,info.xc);
%
if isfield(info,'lc')
    channel_length  = cat(2,channel_length,info.lc);
    channel_amplitude_ratio = cat(2,channel_amplitude_ratio,info.rc);
else
    channel_length(ii)  = 0;
    channel_amplitude_ratio(ii) = 0;
end
%
depFile = [info.rootMat,'funwave_',runDIR,'_dep.nc'];
rotFile = [info.rootMat,'funwave_',runDIR,'_velocity_decomposition.nc'];
momFile = [info.rootMat,'funwave_',runDIR,'_MomentumTerms.nc'];
%
% get (x,y,dep)
x   = ncread(depFile,'x');
y   = ncread(depFile,'y');
dep = ncread(depFile,'dep');
%
% get total exchange velocity:
fileInfo = ncinfo(rotFile);
variableNames = {fileInfo.Variables.Name};
if ~ismember('Urot_mean',variableNames)
    Umean = mean(ncread(momFile,'umean'),3);
    Vmean = mean(ncread(momFile,'vmean'),3);
    [~,Umean,Vmean,~,~,~]=get_vel_decomposition_reGRID(Umean,Vmean,info.dx,info.dy);
else
    Umean = ncread(rotFile,'Urot_mean');
    Vmean = ncread(rotFile,'Vrot_mean');
end
U     = ncread(rotFile,'Urot');
% $$$ V     = ncread(rotFile,'Vrot');
disp('using time-averaged waterlevel... bug in source code')
ETA   = mean(ncread(momFile,'etamean'),3);
% $$$ ETA   = ncread(rotFile,'eta');
% $$$ ETA(ETA>dep) = 0;
%
% create depth mask (min-depth-resolved=0.01m, min-depth-normalize=0.1m)
H         = dep+ETA;
mask      = H>0.01;
H(~mask)  = 0;
Hmean     = max( mean(H,3), 0.1);
%
iOFF      = (Umean+U)>0;
iOFF_mean = (Umean)>0;
iOFF_eddy = (U)>0;
%
VORT_mean = curl(x,y,Umean,Vmean);
VORT      = ncread(rotFile,'VORT');
%
if ii>1
    nx = length(x);
    ne = size(ENS,1);
    if nx>ne
        Uex     (ne:nx,:) = 0;
        Uex_mean(ne:nx,:) = 0;
        Uex_eddy(ne:nx,:) = 0;        
        ENS     (ne:nx,:) = 0;
        ENS_mean(ne:nx,:) = 0;
        ENS_eddy(ne:nx,:) = 0;
        dETA    (ne:nx,:) = 0;
    elseif ne>nx
        U        (:,nx:ne,:) = 0;
        Umean    (:,nx:ne,:) = 0;
        Vmean    (:,nx:ne,:) = 0;
        iOFF     (:,nx:ne,:) = 0;
        iOFF_mean(:,nx:ne,:) = 0;
        iOFF_eddy(:,nx:ne,:) = 0;
        mask     (:,nx:ne,:) = 0;
        H        (:,nx:ne,:) = 0;
        Hmean    (:,nx:ne,:) = 0;
        VORT_mean(:,nx:ne,:) = 0;
        VORT     (:,nx:ne,:) = 0;
        ETA      (:,nx:ne,:) = 0;
    end
end
%
Uex(:,ii)       = mean( (Umean+U).*iOFF     .*H, [1 3]) ./ mean(Hmean, [1 3]);
Uex_mean(:,ii)  = mean( (Umean  ).*iOFF_mean.*H, [1 3]) ./ mean(Hmean, [1 3]);
Uex_eddy(:,ii)  = mean( (U      ).*iOFF_eddy.*H, [1 3]) ./ mean(Hmean, [1 3]);
%
ENS(:,ii)       = mean( (VORT_mean+VORT).^2.*H, [1 3])./mean( Hmean, 1);
ENS_mean(:,ii)  = mean( (VORT_mean     ).^2.*H, [1 3])./mean( Hmean, 1);
ENS_eddy(:,ii)  = mean( (VORT          ).^2.*H, [1 3])./mean( Hmean, 1);
%
%% not all runs have a channel, but for those that do it's located at y~Ly/2
if isfield(info,'lc')
    iY0 =  y<info.Ly/2 - 5*info.lc | y>info.Ly/2 + 5*info.lc;
    iY1 = (y>info.Ly/2 -   info.lc & y<info.Ly/2 +   info.lc);
    iY2 = (y>info.Ly/2 - 3*info.lc & y<info.Ly/2 + 3*info.lc);    
    iX0 = (x>0.5*info.xc & x<info.xc)';
    iX1 = find(x>=info.xc,1,'first');
    %
    Uex_channel     (:,ii)  = mean( (Umean(iY1,:,:)+U(iY1,:,:)).*iOFF(iY1,:,:)     .*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
    Uex_mean_channel(:,ii)  = mean( (Umean(iY1,:,:)           ).*iOFF_mean(iY1,:,:).*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
    Uex_eddy_channel(:,ii)  = mean( (U(iY1,:,:)               ).*iOFF_eddy(iY1,:,:).*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
    %
    Uex_ambient     (:,ii)  = mean( (Umean(iY0,:,:)+U(iY0,:,:)).*iOFF(iY0,:,:)     .*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
    Uex_mean_ambient(:,ii)  = mean( (Umean(iY0,:,:)           ).*iOFF_mean(iY0,:,:).*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
    Uex_eddy_ambient(:,ii)  = mean( (U(iY0,:,:)               ).*iOFF_eddy(iY0,:,:).*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
    %
    ENS_channel     (:,ii)  = mean( (VORT_mean(iY1,:,:)+VORT(iY1,:,:)).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
    ENS_mean_channel(:,ii)  = mean( (VORT_mean(iY1,:,:)              ).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
    ENS_eddy_channel(:,ii)  = mean( (VORT(iY1,:,:)                   ).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
    %
    ENS_ambient     (:,ii)  = mean( (VORT_mean(iY0,:,:)+VORT(iY0,:,:)).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
    ENS_mean_ambient(:,ii)  = mean( (VORT_mean(iY0,:,:)              ).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
    ENS_eddy_ambient(:,ii)  = mean( (VORT(iY0,:,:)                   ).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
    %
    deta= mean(ETA(iY1,:),1)-mean(ETA(iY0,:),1);
    dETA     (:,ii)= deta;
    Uscale   (ii)  = real( sqrt(-2*9.8*mean(deta(iX0))));
    Umax     (ii)  = max(U(iY1,iX1,:)+Umean(iY1,iX1,:), [], [1 3]);
    Umax_mean(ii)  = max(             Umean(iY1,iX1,:), [], [1 3]);
    Umax_eddy(ii)  = max(U(iY1,iX1,:)                 , [], [1 3]);        
    %
    figure(fig0)
    subplot(N,1,ii)
    imagesc(x,(y(iY2)-info.Ly/2)/info.lc,Umean(iY2,:)), caxis([-0.5 0.5]), colormap(cmocean('balance'))
    xline(x(iX1),'--g')
    yline([-1 1],'--b')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[50 400])
    str = split(cbttl,'[');
    annotation('textbox','units','normalized','position',[pos(1:2)+[0 0.25].*pos(3:4), 0.5, 0.1],...
               'string',[str{1},'$=~',cblbl{ii},'$~[',str{2}],...
               'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',6,'backgroundcolor','none')    
    if ii==1
        title(sprintf('max$(\\langle u\\rangle)$: %s',NAME),'interpreter','latex')
        set(gca,'xticklabel',[])
    elseif ii<N
        set(gca,'xticklabel',[])
        if ii==floor(N/2)
            ylabel('$(y-y_0)/L_c$ [m]','interpreter','latex')
        end
    end
    colorbar
    %
    %% make subplots of each RC transport/max estimate
    figure(fig00)
    subplot(N,1,ii)
    imagesc(x,(y(iY2)-info.Ly/2)/info.lc,max(U(iY2,:,:),[],3)), caxis([0 0.5]), colormap(cmocean('amp'))
    xline(x(iX1),'--g')
    yline([-1 1],'--b')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[50 400])
    str = split(cbttl,'[');    
    annotation('textbox','units','normalized','position',[pos(1:2)+[0 0.25].*pos(3:4), 0.5, 0.1],...
               'string',[str{1},'$=~',cblbl{ii},'$~[',str{2}],...
               'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',6,'backgroundcolor','none')    
    if ii==1
        title(sprintf('max$(\\bar{u})$: %s',NAME),'interpreter','latex')
        set(gca,'xticklabel',[])
    elseif ii<N
        set(gca,'xticklabel',[])
        if ii==floor(N/2)
            ylabel('$(y-y_0)/L_c$ [m]','interpreter','latex')
        end
    end
    colorbar
else
    dETA     (1:size(Uex,1),ii)=0;
    Uscale   (ii) = 0;
    Umax     (ii) = 0;
    Umax_mean(ii) = 0;
    Umax_eddy(ii) = 0;
    Uex_channel     (:,ii)  = 0*x;
    Uex_mean_channel(:,ii)  = 0*x;
    Uex_eddy_channel(:,ii)  = 0*x;
    %
    Uex_ambient     (:,ii)  = 0*x;
    Uex_mean_ambient(:,ii)  = 0*x;
    Uex_eddy_ambient(:,ii)  = 0*x;
    %
    ENS_channel     (:,ii)  = 0*x;
    ENS_mean_channel(:,ii)  = 0*x;
    ENS_eddy_channel(:,ii)  = 0*x;
    %
    ENS_ambient     (:,ii)  = 0*x;
    ENS_mean_ambient(:,ii)  = 0*x;
    ENS_eddy_ambient(:,ii)  = 0*x;
end
end

figure(fig0);
xlabel('$x$ [m]','interpreter','latex')
figname = [figDIR,'maximum_channel_mean_speed',NAME,'.pdf'];
exportgraphics(fig0,figname)

figure(fig00);
xlabel('$x$ [m]','interpreter','latex')
figname = [figDIR,'maximum_channel_eddy_speed',NAME,'.pdf'];
exportgraphics(fig00,figname)


cm = cmocean('thermal',N+1);
cm = cm(1:N,:);

% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.5;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
ppos3 = [xm       ym+2*(ph+ag)   pw ph];
cbpos = [xm+pw+ag ym       ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ag+3*ph];

%% Exchange Velocity
fig1 = figure('units','centimeters');
fig1.Position(3:4)=ps;
set(fig1,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,Uex,'-');
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total Exchange Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,Uex_eddy,'-');
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy Exchange Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,Uex_mean,'-');
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean Exchange Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Uex_',NAME,'.pdf'];
exportgraphics(fig1,figname)
% close(fig1)
%%




%% Enstrophy
fig2 = figure('units','centimeters');
fig2.Position(3:4)=ps;
set(fig2,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,sqrt(ENS),'-');
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,sqrt(ENS_eddy),'-');
ylabel('$\langle\omega^2\rangle_{(y,t)}^{1/2}$ [s$^{-1}$]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,sqrt(ENS_mean),'-');
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'ENS_',NAME,'.pdf'];
exportgraphics(fig2,figname)
%%

if sum(dETA~=0,'all')>0
%% Sea-surface gradient
ps1  = [2*xm+pw+6*ag  ym+ag+ph];
fig3 = figure('units','centimeters');
fig3.Position(3:4)=ps1;
set(fig3,'papersize',ps1,'paperposition',[0 0 ps1]);
colororder(cm)

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,dETA,'-');
hold on,xline( [x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))],'--b')
ylabel('$\Delta \eta_y$ [m]','interpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'a) Mean Sealevel Difference:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'dETA_',NAME,'.pdf'];
exportgraphics(fig3,figname)
%%

%% Exchange Velocity
fig4 = figure('units','centimeters');
fig4.Position(3:4)=ps;
set(fig4,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,Uex_channel,'-',x,Uex_ambient,':');
hold on, xline(info.xc,'--g')
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total Exchange Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,Uex_eddy_channel,'-',x,Uex_eddy_ambient,':');
hold on, xline(info.xc,'--g')
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy Exchange Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,Uex_mean_channel,'-',x,Uex_mean_ambient,':');
hold on, xline(info.xc,'--g')
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean Exchange Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Uex_',NAME,'.pdf'];
exportgraphics(fig4,figname)
%%

%% Velocity magnitudes
ppos_eq = [xm ym pw pw];
ps1  = [2*xm+pw+6*ag  ym+ag+pw];
fig5 = figure('units','centimeters');
fig5.Position(3:4)=ps1;
set(fig5,'papersize',ps1,'paperposition',[0 0 ps1]);
colororder(cm)

a1 = axes('units','centimeters','position',ppos_eq);
for ii=1:N
    plot(Uscale(ii),Umax(ii),'o','markerfacecolor',cm(ii,:),'markeredgecolor',cm(ii,:)); hold on
    plot(Uscale(ii),Umax_mean(ii),'s','markeredgecolor',cm(ii,:));
    plot(Uscale(ii),Umax_eddy(ii),'d','markeredgecolor',cm(ii,:));
end
lims = max(xlim,ylim);
lims(1)=0;
xlim(a1,lims), ylim(lims)
hold on,plot(lims,lims,'--k')
ylabel('max$(u)$ [m/s]','interpreter','latex')
xlabel('$\sqrt{-2g\Delta \eta}$ [m/s]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos_eq(1:2)+[0 0.9].*ppos_eq(3:4), 0.3, 0.3],...
           'string',{'a) Rip-channel Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
legend(a1.Children([end,end-1,end-2]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
grid([a1],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Umax_vs_Uscale_',NAME,'.pdf'];
exportgraphics(fig5,figname)
%%

%% Exchange Velocity magnitudes
ps2  = [2*xm+pw+6*ag  ym+2*(ag+ph)];
fig6 = figure('units','centimeters');
fig6.Position(3:4)=ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
% colororder(cm)

a2   = axes('units','centimeters','position',ppos2);
vals = str2num(char(cblbl'));
plot(vals,Uex_channel(iX1,:),'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Uex_mean_channel(iX1,:),'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Uex_eddy_channel(iX1,:),'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
set(a2,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xticklabel',[])
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'a) Rip-channel Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
legend(a2.Children([3,2,1]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
grid([a2],'on')

a1   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Uex_ambient(iX1,:),'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Uex_mean_ambient(iX1,:),'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Uex_eddy_ambient(iX1,:),'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
ylabel('~~~~~~~~~~~~~~~~~~~~~~~~$U_\mathrm{ex}$ [m/s]','interpreter','latex')
xlabel(cbttl,'interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'b) Ambient Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
grid([a1],'on')

figname = [figDIR,filesep,'Uex_channel_vs_ambient_',NAME,'.pdf'];
exportgraphics(fig6,figname)
%%

%% Exchange Velocity magnitudes
fig7 = figure('units','centimeters');
fig7.Position(3:4)=ps1;
set(fig7,'papersize',ps1,'paperposition',[0 0 ps1]);
% colororder(cm)

a1   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Uscale,'x','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umax,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umax_mean,'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Umax_eddy,'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylabel('[m/s]','interpreter','latex')
xlabel(cbttl,'interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
ylims = ylim;
ylim([0 1.2*ylims(2)]);
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'b) Rip-Channel Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
legend(a1.Children([4,3,2,1]),{'$\sqrt{-2g\Delta \eta}$','$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
grid([a1],'on')

figname = [figDIR,filesep,'Umax_',NAME,'.pdf'];
exportgraphics(fig7,figname)
%%

% $$$ legend(a1.Children([1,2,3]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)

% $$$ %% Exchange Velocity magnitudes
% $$$ ppos_eq = [xm ym pw pw];
% $$$ ps1  = [2*xm+pw+6*ag  ym+ag+pw];
% $$$ fig6 = figure('units','centimeters');
% $$$ fig6.Position(3:4)=ps1;
% $$$ set(fig6,'papersize',ps1,'paperposition',[0 0 ps1]);
% $$$ colororder(cm)
% $$$ 
% $$$ a1 = axes('units','centimeters','position',ppos_eq);
% $$$ for ii=1:N
% $$$     tmp1(ii) = Uex_channel(iX1,ii);
% $$$     tmp2(ii) = Uex_mean_channel(iX1,ii);
% $$$     tmp3(ii) = Uex_eddy_channel(iX1,ii);
% $$$     plot(Uscale(ii),tmp1(ii),'o','markerfacecolor',cm(ii,:),'markeredgecolor',cm(ii,:)); hold on
% $$$     plot(Uscale(ii),tmp2(ii),'s','markeredgecolor',cm(ii,:));
% $$$     plot(Uscale(ii),tmp3(ii),'d','markeredgecolor',cm(ii,:));
% $$$ end
% $$$ xlims = xlim;
% $$$ ylims = ylim;
% $$$ xlims(1)=0; ylims(1)=0;
% $$$ xlim(a1,xlims), ylim(a1,ylims)
% $$$ hold on,plot(xlims,xlims/12,'--k')
% $$$ ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
% $$$ xlabel('$\sqrt{-2g\Delta \eta}$ [m/s]','interpreter','latex')
% $$$ set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
% $$$ title(sprintf('%s',NAME),'interpreter','latex')
% $$$ annotation('textbox','units','centimeters','position',[ppos_eq(1:2)+[0 0.9].*ppos_eq(3:4), 0.3, 0.3],...
% $$$            'string',{'a) Rip-channel Velocity Scales: (--) 1/12-slope'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')
% $$$ legend(a1.Children([end,end-1,end-2]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
% $$$ grid([a1],'on')
% $$$ 
% $$$ cb = axes('units','centimeters','position',cbpos);
% $$$ imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
% $$$ xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')
% $$$ 
% $$$ figname = [figDIR,filesep,'Uex_vs_Uscale_channel_',NAME,'.pdf'];
% $$$ exportgraphics(fig6,figname)
% $$$ %%
% $$$ 
% $$$ %% Exchange Velocity magnitudes
% $$$ ppos_eq = [xm ym pw pw];
% $$$ ps1  = [2*xm+pw+6*ag  ym+ag+pw];
% $$$ fig7 = figure('units','centimeters');
% $$$ fig7.Position(3:4)=ps1;
% $$$ set(fig7,'papersize',ps1,'paperposition',[0 0 ps1]);
% $$$ colororder(cm)
% $$$ 
% $$$ a1 = axes('units','centimeters','position',ppos_eq);
% $$$ for ii=1:N
% $$$     tmp1(ii) = Uex_ambient(iX1,ii);
% $$$     tmp2(ii) = Uex_mean_ambient(iX1,ii);
% $$$     tmp3(ii) = Uex_eddy_ambient(iX1,ii);
% $$$     plot(Uscale(ii),tmp1(ii),'o','markerfacecolor',cm(ii,:),'markeredgecolor',cm(ii,:)); hold on
% $$$     plot(Uscale(ii),tmp2(ii),'s','markeredgecolor',cm(ii,:));
% $$$     plot(Uscale(ii),tmp3(ii),'d','markeredgecolor',cm(ii,:));
% $$$ end
% $$$ xlims = xlim;
% $$$ ylims = ylim;
% $$$ xlims(1)=0; ylims(1)=0;
% $$$ xlim(a1,xlims), ylim(a1,ylims)
% $$$ hold on,plot(xlims,xlims/4,'--k')
% $$$ ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
% $$$ xlabel('$\sqrt{-2g\Delta \eta}$ [m/s]','interpreter','latex')
% $$$ set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
% $$$ title(sprintf('%s',NAME),'interpreter','latex')
% $$$ annotation('textbox','units','centimeters','position',[ppos_eq(1:2)+[0 0.9].*ppos_eq(3:4), 0.3, 0.3],...
% $$$            'string',{'a) Ambient Velocity Scales: (--) 1/4-slope'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')
% $$$ legend(a1.Children([end,end-1,end-2]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
% $$$ grid([a1],'on')
% $$$ 
% $$$ cb = axes('units','centimeters','position',cbpos);
% $$$ imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
% $$$ xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')
% $$$ 
% $$$ figname = [figDIR,filesep,'Uex_vs_Uscale_ambient_',NAME,'.pdf'];
% $$$ exportgraphics(fig7,figname)
% $$$ %%
else
    %% Exchange Velocity magnitudes
ps1  = [2*xm+pw+6*ag  ym+1*(ag+ph)];
fig6 = figure('units','centimeters');
fig6.Position(3:4)=ps1;
set(fig6,'papersize',ps1,'paperposition',[0 0 ps1]);
% colororder(cm)

a2   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Uex(iX1,:),'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Uex_mean(iX1,:),'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Uex_eddy(iX1,:),'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
set(a2,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
xlabel(cbttl,'interpreter','latex')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'a) Bar-Crest Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
legend(a2.Children([3,2,1]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
grid([a2],'on')

figname = [figDIR,filesep,'Uex_bar_crest_',NAME,'.pdf'];
exportgraphics(fig6,figname)

end

save([outDIR,'BulkVelocityStats_',NAME,'.mat'],'x','Uex','Uex_mean','Uex_eddy','ENS','ENS_mean','ENS_eddy','dETA','Uscale','Umax','Umax_mean','Umax_eddy','runIDs','run_dirs','cblbl','cbttl','Uex_channel','Uex_ambient','Uex_mean_channel','Uex_mean_ambient','Uex_eddy_channel','Uex_eddy_ambient','ENS_channel','ENS_ambient','ENS_mean_channel','ENS_mean_ambient','ENS_eddy_channel','ENS_eddy_ambient','height','period','spread','direction','bar_width','channel_length','bar_amplitude','channel_amplitude_ratio')
end