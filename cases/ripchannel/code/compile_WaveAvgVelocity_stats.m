%% create compendium plot of exchange and enstrophy for a set of runs:
clear all
close all
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'

%% 1) need a list of run-directories: 'uniRip-ter2D','uniRip-bar2D','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10',
NAMES    = {'spreadRip-barRip1'};%{'spreadRip-barRip0','spreadRip-barRip1','spreadRip-terRip1','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10'};

for nn = 1:length(NAMES)
    clearvars -except rootDIR figDIR outDIR NAMES nn
    
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
N     = length(runIDs);
%
% for velocity plots (maximum of 5-panels)
xm = 2;
ym = 2;
pw = 2.5;% 400 in x
ph = 6;% +/- 500 in y
ag = 0.5;
ppos1 = [xm ym pw ph];
ppos2 = [xm+pw+ag ym pw ph];
ppos3 = [xm+2*(pw+ag) ym pw ph];
ppos4 = [xm+3*(pw+ag) ym pw ph];
ppos5 = [xm+4*(pw+ag) ym pw ph];
ps    = [2*xm+N*(pw+ag)+6*ag 2*ym+ph];
cbpos = [xm+N*(pw+ag)+ag, ym, ag, ph/2];
%
fig0  = figure('units','centimeters');
fig00 = figure('units','centimeters');
fig000 = figure('units','centimeters');
fig0000 = figure('units','centimeters');
fig00000 = figure('units','centimeters');
fig0.Position(3:4)  = ps;
fig00.Position(3:4) = ps;
fig000.Position(3:4) = ps;
fig0000.Position(3:4) = ps;
fig00000.Position(3:4) = ps;
set(fig0 ,'papersize',ps,'paperposition',[0 0 ps])
set(fig00,'papersize',ps,'paperposition',[0 0 ps])
set(fig000,'papersize',ps,'paperposition',[0 0 ps])
set(fig0000,'papersize',ps,'paperposition',[0 0 ps])
set(fig00000,'papersize',ps,'paperposition',[0 0 ps])
%
% $$$ fig1  = figure('units','centimeters');
% $$$ fig11 = figure('units','centimeters');
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
h0  = dep(1,:);
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
V     = ncread(rotFile,'Vrot');
%
t     = ncread(rotFile,'t');
% $$$ V     = ncread(rotFile,'Vrot');
disp('using time-averaged waterlevel... bug in source code')
ETAmean   = mean(ncread(momFile,'etamean'),3);
ETA       = ETAmean;
% $$$ ETA   = ncread(rotFile,'eta');
% $$$ ETA(ETA>dep) = 0;
%
%% create depth mask (min-depth-resolved=0.01m, min-depth-normalize=0.1m)
H         = dep+ETA;
mask      = H>0.01;
H(~mask)  = 0;
Hmean     = max( mean(H,3), 0.1);
%
%% logical array fro offshore flow
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
        Uke     (ne:nx,:) = 0;
        Ueke    (ne:nx,:) = 0;
        Umke    (ne:nx,:) = 0;
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
    %
    %
    if exist('Umax_vs_t','var')
        nt1 = size(Umax_vs_t,1);
        nt2 = size(U,3);
        if nt1<nt2
            Umax_vs_t(nt1:nt2,:)=nan;
            Umax_eddy_vs_t(nt1:nt2,:)=nan;            
        elseif nt2<nt1
            U(:,:,nt2:nt1) = nan;
            if size(H,3)>1
                H(:,:,nt2:nt1) = nan;
            end
            VORT     (:,:,nt2:nt1) = nan;
            iOFF     (:,:,nt2:nt1) = 0;            
            iOFF_eddy(:,:,nt2:nt1) = 0;
        end
    end
end
%
%% Cross-shore transports:
T     = (Umean+U).*H; T    (~iOFF) = nan;
Tmean = Umean.*Hmean; Tmean(~iOFF_mean) = nan;
Teddy = U.*H;         Teddy(~iOFF_eddy) = nan;
%
%% Full alongshore domain statistics
Uex(:,ii)       = mean( sum(T    ,1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1 ,'omitnan');
Uex_mean(:,ii)  = mean( sum(Tmean,1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1 ,'omitnan');
Uex_eddy(:,ii)  = mean( sum(Teddy,1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1 ,'omitnan');
%
Uke(:,ii)   = sqrt(mean(sum( (Umean+U).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
Umke(:,ii)  = sqrt(mean(sum( (Umean  ).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
Ueke(:,ii)  = sqrt(mean(sum( (U      ).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
%
ENS(:,ii)       = mean( (VORT_mean+VORT).^2.*H, [1 3],'omitnan')./mean( Hmean, [1 3],'omitnan');
ENS_mean(:,ii)  = mean( (VORT_mean     ).^2.*H, [1 3],'omitnan')./mean( Hmean, [1 3],'omitnan');
ENS_eddy(:,ii)  = mean( (VORT          ).^2.*H, [1 3],'omitnan')./mean( Hmean, [1 3],'omitnan');
%
%
iX0 = (x>0.6*info.xc & x<0.9*info.xc)';
iX1 = find(x>=info.xc,1,'first');
iX2 = find(x>=info.xc-10 & x<=info.xc+10);
iX3 = find(x>=2*info.xc-10,1,'first');
%
%
%% not all runs have a channel, but for those that do it's located at y~Ly/2
if isfield(info,'lc')
    iY0 =  y<info.Ly/2 - 5*info.lc | y>info.Ly/2 + 5*info.lc;
    iY1 = (y>info.Ly/2 - 1*info.lc & y<info.Ly/2 + 1*info.lc);
    iY2 = (y>info.Ly/2 - 2*info.lc & y<info.Ly/2 + 2*info.lc);    
    %
    Uex_channel     (:,ii)  = mean( sum(T    (iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:,:), 1,'omitnan');
    Uex_mean_channel(:,ii)  = mean( sum(Tmean(iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:,:), 1,'omitnan');
    Uex_eddy_channel(:,ii)  = mean( sum(Teddy(iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:,:), 1,'omitnan');
    %
    Uex_ambient     (:,ii)  = mean( sum(T    (iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:,:), 1,'omitnan');
    Uex_mean_ambient(:,ii)  = mean( sum(Tmean(iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:,:), 1,'omitnan');
    Uex_eddy_ambient(:,ii)  = mean( sum(Teddy(iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:,:), 1,'omitnan');
    %
    ENS_channel     (:,ii)  = mean( (VORT_mean(iY2,:,:)+VORT(iY2,:,:)).^2.*H(iY2,:,:), [1 3],'omitnan')./mean( Hmean(iY2,:,:), [1 3],'omitnan');
    ENS_mean_channel(:,ii)  = mean( (VORT_mean(iY2,:,:)              ).^2.*H(iY2,:,:), [1 3],'omitnan')./mean( Hmean(iY2,:,:), [1 3],'omitnan');
    ENS_eddy_channel(:,ii)  = mean( (VORT(iY2,:,:)                   ).^2.*H(iY2,:,:), [1 3],'omitnan')./mean( Hmean(iY2,:,:), [1 3],'omitnan');
    %
    ENS_ambient     (:,ii)  = mean( (VORT_mean(iY0,:,:)+VORT(iY0,:,:)).^2.*H(iY0,:,:), [1 3],'omitnan')./mean( Hmean(iY0,:,:), [1 3],'omitnan');
    ENS_mean_ambient(:,ii)  = mean( (VORT_mean(iY0,:,:)              ).^2.*H(iY0,:,:), [1 3],'omitnan')./mean( Hmean(iY0,:,:), [1 3],'omitnan');
    ENS_eddy_ambient(:,ii)  = mean( (VORT(iY0,:,:)                   ).^2.*H(iY0,:,:), [1 3],'omitnan')./mean( Hmean(iY0,:,:), [1 3],'omitnan');
    %
    Uke_channel (:,ii)  = sqrt(mean(sum( (Umean(iY2,:,:)+U(iY2,:,:)).^2.*H(iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:), 1,'omitnan'));
    Umke_channel(:,ii)  = sqrt(mean(sum( (Umean(iY2,:,:)           ).^2.*H(iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:), 1,'omitnan'));
    Ueke_channel(:,ii)  = sqrt(mean(sum( (U(iY2,:,:)               ).^2.*H(iY2,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY2,:), 1,'omitnan'));
    %
    Uke_ambient (:,ii)  = sqrt(mean(sum( (Umean(iY0,:,:)+U(iY0,:,:)).^2.*H(iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:), 1,'omitnan'));
    Umke_ambient(:,ii)  = sqrt(mean(sum( (Umean(iY0,:,:)           ).^2.*H(iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:), 1,'omitnan'));
    Ueke_ambient(:,ii)  = sqrt(mean(sum( (U(iY0,:,:)               ).^2.*H(iY0,:,:),1,'omitnan'), 3,'omitnan') ./ sum(Hmean(iY0,:), 1,'omitnan'));
    %
    %% Potential difference:
    deta= mean(ETAmean(iY1,:),1)-mean(ETAmean(iY0,:),1);
    dETA     (:,ii)= deta;
    Uscale   (ii)  = real( sqrt(-2*9.8*mean(deta(iX0))));
    %
    %% Alongshore velocity stats:
    ETA_vs_y    (:,ii) = mean(ETAmean(:,iX0,:), [2 3],'omitnan');
    Vmean_feeder(:,ii) = mean(Vmean(:,iX0,:).*Hmean(:,iX0,:), [2 3])./mean(Hmean(:,iX0,:), [2 3]);
    Veddy_feeder(:,ii) = sqrt( sum (V(:,iX0,:).^2.*H(:,iX0,:), [2 3])./sum(H(:,iX0,:), [2 3]));
    iYplus     = find(y>info.Ly/2 & y<info.Ly/2+6*info.lc);
    iYminus    = find(y<info.Ly/2 & y>info.Ly/2-6*info.lc);    
    Vscale(ii) = 0.5*( -min(Vmean_feeder(iYplus,ii),[],1) + max(Vmean_feeder(iYminus,ii),[],1) );
    Vrms  (ii) = rms(Vmean_feeder(iYplus | iYminus,ii),1);
    %
    %% Use cross-shore kinetic energy to locate the outer-surfzone "maximum":
    iOuterSZ = (x>0.9*info.xc & x<1.5*info.xc);
    [Uke_channel_max(ii) , idx_ke_max ] = max(Uke_channel (:,ii).*iOuterSZ);
    [Umke_channel_max(ii), idx_mke_max] = max(Umke_channel(:,ii).*iOuterSZ);
    [Ueke_channel_max(ii), idx_eke_max] = max(Ueke_channel(:,ii).*iOuterSZ);
    Xke_max (ii)  = x(idx_ke_max);
    Xmke_max(ii) = x(idx_mke_max);
    Xeke_max(ii) = x(idx_eke_max);
    %
    % time-series of "maximum" speed at cross-shore ke maximum
    Umax_vs_t      (:,ii) = max(Umean(iY2,idx_ke_max) + U(iY2,idx_ke_max,:),[],1);
    Umax_eddy_vs_t (:,ii) = max(U(iY2,idx_ke_max,:),[],1);        
    %
    % find the time where the maximum velocity occurs in simulation:
    [~,idt_ke_max] = max( max(Umean(iY2,idx_ke_max) + U(iY2,idx_ke_max,:),[],1),[],3,'omitnan');
    Umax_vs_y  (:,ii) = Umean(:,idx_ke_max) + U(:,idx_ke_max,idt_ke_max);
    Umean_vs_y (:,ii) = Umean(:,idx_ke_max);
    Umean_vs_y_offshore(:,ii) = Umean(:,iX3,:);
    Urms_eddy_vs_y  (:,ii) = rms(U(:,idx_ke_max,:),3,'omitnan');
    %
    Umax     (ii)  = mean( max(Umean(iY1,idx_ke_max) + U(iY1,idx_ke_max,:), [], 1),3,'omitnan');
    Umax_mean(ii)  =       max(Umean(iY1,idx_ke_max)   , [], 1 ,'omitnan');
    Umean_eddy(ii) = mean( max(    U(iY1,idx_ke_max,:) , [], 1 ,'omitnan'),3,'omitnan');
    Urms_eddy(ii)  = rms ( max(    U(iY1,idx_ke_max,:) , [], 1 ,'omitnan'),3,'omitnan');
    %
    [Uke_ambient_max(ii) , ~ ] = max(Uke_ambient (:,ii).*iOuterSZ);
    [Umke_ambient_max(ii), ~ ] = max(Umke_ambient(:,ii).*iOuterSZ);
    [Ueke_ambient_max(ii), ~ ] = max(Ueke_ambient(:,ii).*iOuterSZ);
    %
    figure(fig0)
    %subplot(N,1,ii)
    eval(sprintf('pos = ppos%d;',ii))
    axes('units','centimeters','position',pos);
    imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,Umean), caxis([-0.5 0.5]), colormap(cmocean('balance'))
    hold on,
    contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
    xline(Xke_max(ii),'--m')
    yline([-2 2],'--g')
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',0.5*ph/pw*450/info.lc*[-1 1])
    str = split(cbttl,'[');
    title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
    if ii==1
        ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
    else
        set(gca,'yticklabel',[])
        if ii==floor(N/2)
            xlabel('$(x-x_sl)/X_c$ [~]','interpreter','latex')
        end
    end
    %
    %% make subplots of each RC transport/max estimate
    figure(fig00)
    axes('units','centimeters','position',pos);
    imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,rms(U,3,'omitnan')), caxis([0 0.5]), colormap(cmocean('amp'))
    hold on,
    contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
    xline(Xke_max(ii),'--m')
    yline([-2 2],'--g')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',0.5*ph/pw*450/info.lc*[-1 1])
    str = split(cbttl,'[');    
    title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
    if ii==1
        ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
    else
        set(gca,'yticklabel',[])
        if ii==floor(N/2)
            xlabel('$(x-x_sl)/X_c$ [~]','interpreter','latex')
        end
    end
    %
    %%
    figure(fig000)
    eval(sprintf('pos = ppos%d;',ii))
    axes('units','centimeters','position',pos);
    imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,Vmean), caxis([-0.5 0.5]), colormap(cmocean('balance'))
    hold on,
    contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
    xline( [x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))],'--c')
    yline([-1 1],'--b')
    yline([-5 5],'--c')
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',0.5*ph/pw*450/info.lc*[-1 1])
    str = split(cbttl,'[');
    title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
    if ii==1
        ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
    else
        set(gca,'yticklabel',[])
        if ii==floor(N/2)
            xlabel('$(x-x_sl)/X_c$ [~]','interpreter','latex')
        end
    end
    %
    %
    %% make subplots of mean alongshore current
    figure(fig0000)
    axes('units','centimeters','position',pos);
    imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,rms(V,3,'omitnan')), caxis([0 0.5]), colormap(cmocean('amp'))
    hold on,
    contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
    xline( [x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))],'--c')
    yline([-1 1],'--b')
    pos = get(gca,'position');
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',0.5*ph/pw*450/info.lc*[-1 1])
    str = split(cbttl,'[');    
    title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
    if ii==1
        ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
    else
        set(gca,'yticklabel',[])
        if ii==floor(N/2)
            xlabel('$(x-x_sl)/X_c$ [~]','interpreter','latex')
        end
    end
    %
    %
    %% Vorticity plot
    figure(fig00000)
    eval(sprintf('pos = ppos%d;',ii))
    axes('units','centimeters','position',pos);
    imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,VORT_mean), caxis([-0.025 0.025]), colormap(cmocean('curl'))
    hold on,
    [~,h] = contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5);
    %% add quiver
    spanx = round(20/info.dx);
    spany = round(40/info.dy);
    %% want to start from middle of domain:
% $$$     iYquiv = round(info.Ly/2)+ [-40*spany:spany:40*spany]; iYquiv = iYquiv(iYquiv>0 & iYquiv<length(y));
% $$$     [xx,yy] = meshgrid((x(1:spanx:end)-50)/(info.xc-50),(y(iYquiv)-info.Ly/2)/info.lc);
% $$$     uu = Umean(iYquiv,1:spanx:end);
% $$$     vv = Vmean(iYquiv,1:spanx:end);
    [xx,yy] = meshgrid((x(1:spanx:end)-50)/(info.xc-50),(y(1:spany:end)-info.Ly/2)/info.lc);
    uu = Umean(1:spany:end,1:spanx:end);
    vv = Vmean(1:spany:end,1:spanx:end);
    xlims = [0 450]/(info.xc-50);
    ylims = 0.5*(ph/pw)*(400/info.lc)*[-1 1];
    scalex = 30*(info.xc-50)/info.lc/(info.xc-50);
    scaley = 30*(diff(ylims)/diff(xlims))*(pw/ph)/(info.lc);
    quiver(xx(:),yy(:),uu(:)*scalex, vv(:)*scaley,'ShowArrowHead','off','Marker','.','markersize',1)
    quiver(1.8,4,0.5*scalex, 0,'ShowArrowHead','off','Marker','.','markersize',1,'color','k','linewidth',1)
    text(  1.8,4.75,'0.5 m/s','interpreter','latex','fontsize',5)
    %%
    set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
    str = split(cbttl,'[');
    title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
    if ii==1
        ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
    else
        set(gca,'yticklabel',[])
        if ii==floor(N/2)
            xlabel('$(x-x_sl)/X_c$ [~]','interpreter','latex')
        end
    end
    %

else
    dETA     (1:size(Uex,1),ii)=0;
    Uscale   (ii) = 0;
    Umax     (ii) = 0;
    Umax_mean(ii) = 0;
    Umean_eddy(ii) = 0;
    Urms_eddy(ii) = 0;    
    Uex_channel     (:,ii)  = 0*x;
    Uex_mean_channel(:,ii)  = 0*x;
    Uex_eddy_channel(:,ii)  = 0*x;
    %
    Uex_ambient     (:,ii)  = 0*x;
    Uex_mean_ambient(:,ii)  = 0*x;
    Uex_eddy_ambient(:,ii)  = 0*x;
    %
    Uke (:,ii)  = 0*x;
    Umke(:,ii)  = 0*x;
    Ueke(:,ii)  = 0*x;
    %
    Uke_channel (:,ii)  = 0*x;
    Umke_channel(:,ii)  = 0*x;
    Ueke_channel(:,ii)  = 0*x;
    %
    Uke_ambient (:,ii)  = 0*x;
    Umke_ambient(:,ii)  = 0*x;
    Ueke_ambient(:,ii)  = 0*x;
    %
    ENS_channel     (:,ii)  = 0*x;
    ENS_mean_channel(:,ii)  = 0*x;
    ENS_eddy_channel(:,ii)  = 0*x;
    %
    ENS_ambient     (:,ii)  = 0*x;
    ENS_mean_ambient(:,ii)  = 0*x;
    ENS_eddy_ambient(:,ii)  = 0*x;
    %
    %% Use cross-shore kinetic energy to locate the outer-surfzone "maximum":
    Xke_max(ii)  = 0;
    Xmke_max(ii) = 0;
    Xeke_max(ii) = 0;
    %
    Umax_vs_y      (:,ii) = 0*y;
    Umean_vs_y (:,ii) = 0*y;
    Umean_vs_y_offshore (:,ii) = 0*y;    
    Urms_eddy_vs_y (:,ii) = 0*y;
    Xke_max(ii)  = 0;
    Xmke_max(ii) = 0;
    Xeke_max(ii) = 0;
end
end

if isfield(info,'lc')
figure(fig0);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel('$\langle u\rangle$ [m/s]','interpreter','latex')

figname = [figDIR,'mean_x_speed_',NAME,'.pdf'];
exportgraphics(fig0,figname)

figure(fig00);
cm    = cmocean('amp');
clims = [0 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel('rms$(\bar{u})$ [m/s]','interpreter','latex')
% $$$ xlabel('$x$ [m]','interpreter','latex')
figname = [figDIR,'rms_x_speed_',NAME,'.pdf'];
exportgraphics(fig00,figname)
%%
figure(fig000);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel('$\langle v\rangle$ [m/s]','interpreter','latex')

figname = [figDIR,'mean_y_speed_',NAME,'.pdf'];
exportgraphics(fig000,figname)
%
figure(fig0000);
cm    = cmocean('amp');
clims = [0 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel('rms$(\bar{v})$ [m/s]','interpreter','latex')
% $$$ xlabel('$x$ [m]','interpreter','latex')
figname = [figDIR,'rms_y_speed_',NAME,'.pdf'];
exportgraphics(fig0000,figname)
%

%
figure(fig00000);
cm    = cmocean('curl');
clims = [-0.025 0.025];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel('$\langle\omega\rangle$ [1/s]','interpreter','latex')
% $$$ xlabel('$x$ [m]','interpreter','latex')
figname = [figDIR,'vorticity_',NAME,'.pdf'];
exportgraphics(fig00000,figname)
%

end

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

figname = [figDIR,filesep,'Uex_',NAME,'.png'];
exportgraphics(fig1,figname)
% close(fig1)
%%

%% RMS Velocity
fig11 = figure('units','centimeters');
fig11.Position(3:4)=ps;
set(fig11,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,Uke,'-');
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total Cross-shore Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,Ueke,'-');
ylabel('$\mathrm{rms}(u)$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy Cross-shore Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,Umke,'-');
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean Cross-shore Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Urms_',NAME,'.png'];
exportgraphics(fig11,figname)
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

figname = [figDIR,filesep,'ENS_',NAME,'.png'];
exportgraphics(fig2,figname)
%%

if sum(dETA~=0,'all')>0
%% Sea-surface gradient
ps1  = [2*xm+pw+6*ag  ym+ag+ph];
fig3 = figure('units','centimeters');
fig3.Position(3:4)=ps1;
set(fig3,'papersize',ps1,'paperposition',[0 0 ps1]);
colororder(cm)
g = 9.8;
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,g*dETA,'-');
hold on,xline( [x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))],'--c')
ylabel('$g\Delta \eta_y$ [m/s]$^2$','interpreter','latex')
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

figname = [figDIR,filesep,'dETA_',NAME,'.png'];
exportgraphics(fig3,figname)
%%

%% Exchange Velocity
fig4 = figure('units','centimeters');
fig4.Position(3:4)=ps;
set(fig4,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot((x-50)/(info.xc-50),Uex_channel,'-',(x-50)/(info.xc-50),Uex_ambient,':');
% $$$ hold on, xline(Xke_max,'--')
ylims = [0 0.15];
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out','ylim',ylims)
title(sprintf('%s',NAME),'interpreter','latex')
% $$$ annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
% $$$            'string',{'a) Total Exchange Velocity: (-) Channel, (:) Ambient'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')
% $$$ ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot((x-50)/(info.xc-50),Uex_eddy_channel,'-',(x-50)/(info.xc-50),Uex_eddy_ambient,':');
% $$$ hold on, xline(Xke_max,'--')
ylabel('$U_\mathrm{ex}$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
% $$$ annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
% $$$            'string',{'b) Eddy Exchange Velocity: (-) Channel, (:) Ambient'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot((x-50)/(info.xc-50),Uex_mean_channel,'-',(x-50)/(info.xc-50),Uex_mean_ambient,':');
% $$$ hold on, xline(Xke_max,'--')
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
hold on, pLeg = plot(xlim,-999*[1 1],'-k',xlim,-999*[1 1],'--k')
legend(pLeg,'Rip-Channel','Ambient')
% $$$ annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
% $$$            'string',{'c) Mean Exchange Velocity: (-) Channel, (:) Ambient'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Uex_channel_vs_ambient_',NAME,'.pdf'];
exportgraphics(fig4,figname)
%%

%% RMS-Cross-shore Velocity
fig4 = figure('units','centimeters');
fig4.Position(3:4)=ps;
set(fig4,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,Uke_channel,'-',x,Uke_ambient,':');
hold on, xline(Xke_max,'--')
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total Cross-shore Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,Ueke_channel,'-',x,Ueke_ambient,':');
hold on, xline(Xke_max,'--')
ylabel('$\mathrm{rms}(u)$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy Cross-shore Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,Umke_channel,'-',x,Umke_ambient,':');
hold on, xline(Xke_max,'--')
xlabel('$x$ [m]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean Cross-shore Velocity: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Urms_channel_vs_ambient_',NAME,'.png'];
exportgraphics(fig4,figname)
%%

%% Enstrophy
fig2 = figure('units','centimeters');
fig2.Position(3:4)=ps;
set(fig2,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,sqrt(ENS_channel),'-',x,sqrt(ENS_ambient),':');
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Total: (-) Channel, (:) Ambient'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3);

a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x,sqrt(ENS_eddy_channel),'-',x,sqrt(ENS_eddy_ambient),':');
ylabel('$\langle\omega^2\rangle_{(y,t)}^{1/2}$ [s$^{-1}$]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Eddy:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,sqrt(ENS_mean_channel),'-',x,sqrt(ENS_mean_ambient),':');
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

figname = [figDIR,filesep,'ENS_channel_vs_ambient_',NAME,'.png'];
exportgraphics(fig2,figname)
%%


%% Velocity Magnitude versus y direction
fig11 = figure('units','centimeters');
fig11.Position(3:4)=ps;
set(fig11,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot((y-info.Ly/2)/info.lc,Umax_vs_y,'-');
xline([-2 2],'--g')
ylabel('$\mathrm{max}(u(t))$ [m/s]','interpreter','latex')
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out','xlim',[-8 8])
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Max Total Velocity at Peak in rms($u$):'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
ylims = ylim(a3)/2;

a2 = axes('units','centimeters','position',ppos2);
p2 = plot((y-info.Ly/2)./info.lc,Urms_eddy_vs_y,'-');
xline([-2 2],'--g')
ylabel('$\mathrm{rms}(\bar{u}(t))$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'ylim',ylims,'tickdir','out','xlim',[-8 8])
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) RMS-Eddy Cross-shore Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot((y-info.Ly/2)./info.lc,Umean_vs_y,'-');
xline([-2 2],'--g')
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
ylabel('$\langle{u}\rangle$ [m/s]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out','xlim',[-8 8])
% $$$ annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
% $$$            'string',{'c) Mean Cross-shore Velocity:'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Umax_vs_y_',NAME,'.pdf'];
exportgraphics(fig11,figname)


%% Velocity Magnitude versus y direction
ps2  = [2*xm+pw+6*ag  ym+2*(ag+ph)];
fig11 = figure('units','centimeters');
fig11.Position(3:4)=ps2;
set(fig11,'papersize',ps2,'paperposition',[0 0 ps2]);
colororder(cm)

a2 = axes('units','centimeters','position',ppos2);
p2 = plot((y-info.Ly/2)./info.lc,Umean_vs_y_offshore,'-');
xline([-2 2],'--g')
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
ylabel('$\langle{u}\rangle$ [m/s]','interpreter','latex')
xline([-2 2],'--g')
ylabel('$\mathrm{rms}(\bar{u}(t))$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out','xlim',[-8 8])
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) Mean Cross-shore Velocity at 2$x_c$:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
ylims = 1.3*get(a2,'ylim');
a1 = axes('units','centimeters','position',ppos1);
p1 = plot((y-info.Ly/2)./info.lc,Umean_vs_y,'-');
xline([-2 2],'--g')
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
ylabel('$\langle{u}\rangle$ [m/s]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'ylim',ylims,'tickdir','out','xlim',[-8 8])
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'c) Mean Cross-shore Velocity at peak rms$(u)$:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
grid([a1 a2],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Umean_vs_y_offshore_',NAME,'.png'];
exportgraphics(fig11,figname)


%% Alongshore stats versus y direction
fig11 = figure('units','centimeters');
fig11.Position(3:4)=ps;
set(fig11,'papersize',ps,'paperposition',[0 0 ps]);
colororder(cm)

a3 = axes('units','centimeters','position',ppos3);
p3 = plot((y-info.Ly/2)/info.lc,ETA_vs_y*9.8,'-');
xline([-1 1],'--b')
ylabel('$g\langle{\eta}\rangle$ [m/s]$^2$','interpreter','latex')
set(a3,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out','xlim',[-10 10])
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'a) Mid-Surfzone Sealevel:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
% $$$ ylims = ylim(a3)/2;

a2 = axes('units','centimeters','position',ppos2);
p2 = plot((y-info.Ly/2)./info.lc,Veddy_feeder,'-');
ylabel('$\mathrm{rms}(\bar{v}(t))$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xticklabel',[],'fontsize',10,'tickdir','out','xlim',[-10 10])
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'b) RMS-Eddy Alongshore Velocity:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    

a1 = axes('units','centimeters','position',ppos1);
p1 = plot((y-info.Ly/2)./info.lc,Vmean_feeder,'-');
xline([-5 5],'--c')
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
ylabel('$\langle{v}\rangle$ [m/s]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xlim',[-10 10])
% $$$ annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
% $$$            'string',{'c) Mean Alongshore Velocity:'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',8,'backgroundcolor','none')    
grid([a1 a2 a3],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Sealevel_and_feeder_vs_y_',NAME,'.pdf'];
exportgraphics(fig11,figname)



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
    plot(Uscale(ii),Urms_eddy(ii),'d','markeredgecolor',cm(ii,:));
end
lims = max(xlim,ylim);
lims(1)=0;
xlim(a1,lims), ylim(lims)
hold on,plot(lims,lims,'--k')
ylabel('[m/s]','interpreter','latex')
xlabel('$\sqrt{-2g\Delta \eta}$ [m/s]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos_eq(1:2)+[0 0.9].*ppos_eq(3:4), 0.3, 0.3],...
           'string',{'a) Rip-channel Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
legend(a1.Children([end,end-1,end-2]),{'max($\langle{u}\rangle + \bar{u}$)','$\langle{u}\rangle$','rms($\bar{u}$)'},'interpreter','latex','location','southeast','fontsize',8)
grid([a1],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',cblbl,'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,cbttl,'interpreter','latex','horizontalalignment','left')

figname = [figDIR,filesep,'Umax_vs_Uscale_',NAME,'.png'];
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
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'b) Ambient Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
grid([a1],'on')

figname = [figDIR,filesep,'Uex_channel_and_ambient_vs_waves',NAME,'.pdf'];
exportgraphics(fig6,figname)
%%

%% Rip-Velocity magnitudes
fig7 = figure('units','centimeters');
fig7.Position(3:4)=ps1;
set(fig7,'papersize',ps1,'paperposition',[0 0 ps1]);
% colororder(cm)

a1   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Uscale,'x','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umax,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umax_mean,'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Umean_eddy,'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
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

figname = [figDIR,filesep,'Umax_vs_waves_',NAME,'.png'];
exportgraphics(fig7,figname)
%%


%% Alongshore Velocity magnitudes
ps2  = [2*xm+pw+6*ag  ym+2*(ag+ph)];
fig6 = figure('units','centimeters');
fig6.Position(3:4)=ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
% colororder(cm)

a2   = axes('units','centimeters','position',ppos2);
vals = str2num(char(cblbl'));
plot(vals,Uscale,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
ylabel('$\sqrt{-2g\Delta \eta}$ [m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xticklabel',[])
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'a) Sealevel-Based PE-Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')

a1   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Vscale,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
ylabel('$V$ [m/s]','interpreter','latex')
xlabel(cbttl,'interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'b) Alongshore Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
grid([a1],'on')

figname = [figDIR,filesep,'Uscale_and_Vscale_vs_waves',NAME,'.png'];
exportgraphics(fig6,figname)
%%


%% Rip-channel RMS-Velocity magnitudes
ps2  = [2*xm+pw+6*ag  ym+2*(ag+ph)];
fig6 = figure('units','centimeters');
fig6.Position(3:4)=ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
% colororder(cm)

a2   = axes('units','centimeters','position',ppos2);
vals = str2num(char(cblbl'));
plot(vals,Uke_channel_max,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umke_channel_max,'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Ueke_channel_max,'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
ylabel('[m/s]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xticklabel',[])
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'a) Rip-Channel RMS-Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')

a1   = axes('units','centimeters','position',ppos1);
vals = str2num(char(cblbl'));
plot(vals,Uke_ambient_max,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(vals,Umke_ambient_max,'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(vals,Ueke_ambient_max,'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on
ylims = ylim;
ylim([0 1.2*ylims(2)]);
ylabel('[m/s]','interpreter','latex')
xlabel(cbttl,'interpreter','latex')
set(a1,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out')
title(sprintf('%s',NAME),'interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'b) Ambient RMS-Velocity Scales:'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')
grid([a1],'on')
legend({'$\langle{u}\rangle+\bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)

figname = [figDIR,filesep,'Urms_vs_waves',NAME,'.png'];
exportgraphics(fig6,figname)
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
% $$$ figname = [figDIR,filesep,'Uex_vs_Uscale_channel_',NAME,'.png'];
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
% $$$ figname = [figDIR,filesep,'Uex_vs_Uscale_ambient_',NAME,'.png'];
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

figname = [figDIR,filesep,'Uex_bar_crest_',NAME,'.png'];
exportgraphics(fig6,figname)
%%
end

close all
save([outDIR,'BulkVelocityStats_',NAME,'.mat'],'-v7.3','x','y','Uex','Uex_mean','Uex_eddy','ENS','ENS_mean','ENS_eddy','dETA','Uscale','Umax','Umax_mean','Umean_eddy','runIDs','run_dirs','cblbl','cbttl','Uex_channel','Uex_ambient','Uex_mean_channel','Uex_mean_ambient','Uex_eddy_channel','Uex_eddy_ambient','ENS_channel','ENS_ambient','ENS_mean_channel','ENS_mean_ambient','ENS_eddy_channel','ENS_eddy_ambient','height','period','spread','direction','bar_width','channel_length','bar_location','bar_amplitude','channel_amplitude_ratio','Uke','Ueke','Umke','Uke_channel','Ueke_channel','Umke_channel','Uke_ambient','Ueke_ambient','Umke_ambient','Umax_vs_y','Umean_vs_y','Urms_eddy_vs_y','Umax_vs_t','Umax_eddy_vs_t','ETA_vs_y','Vmean_feeder','Veddy_feeder','Uke_channel_max','Umke_channel_max','Ueke_channel_max','Uke_ambient_max','Umke_ambient_max','Ueke_ambient_max','Xke_max','Xmke_max','h0','t')

end