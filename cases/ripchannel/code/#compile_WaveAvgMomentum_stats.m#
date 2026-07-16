%% create compendium plot of momentum terms for a set of runs:
clear all
close all
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'

%% 1) need a list of run-directories: 'uniRip-ter2D','uniRip-bar2D','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10',
NAMES    = {'spreadRip-barRip0','spreadRip-barRip1','spreadRip-terRip1','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10'};

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
%% need figures for:
%      X-Momentum: -PgrdX/H, -DxSxx/H-DySxy/H+BrkDissX/H, -FrcX/H, DxUU+DyUV, Dxuu+Dyuv,
%      Y-Momentum: -PgrdY/H, -DySyy/H-DxSxy/H+BrkDissY/H, -FrcX/H, DyVV+DxUV, Dyvv+Dxuv,
figPGX  = figure('units','centimeters');
figRSX = figure('units','centimeters');
figFRCX = figure('units','centimeters');
figADXmean = figure('units','centimeters');
figADXeddy = figure('units','centimeters');
figPGX.Position(3:4)  = ps;
figRSX.Position(3:4) = ps;
figFRCX.Position(3:4) = ps;
figADXmean.Position(3:4) = ps;
figADXeddy.Position(3:4) = ps;
set(figPGX ,'papersize',ps,'paperposition',[0 0 ps])
set(figRSX,'papersize',ps,'paperposition',[0 0 ps])
set(figFRCX,'papersize',ps,'paperposition',[0 0 ps])
set(figADXmean,'papersize',ps,'paperposition',[0 0 ps])
set(figADXeddy,'papersize',ps,'paperposition',[0 0 ps])
%
figPGY  = figure('units','centimeters');
figRSY = figure('units','centimeters');
figFRCY = figure('units','centimeters');
figADYmean = figure('units','centimeters');
figADYeddy = figure('units','centimeters');
figPGY.Position(3:4)  = ps;
figRSY.Position(3:4) = ps;
figFRCY.Position(3:4) = ps;
figADYmean.Position(3:4) = ps;
figADYeddy.Position(3:4) = ps;
set(figPGY ,'papersize',ps,'paperposition',[0 0 ps])
set(figRSY,'papersize',ps,'paperposition',[0 0 ps])
set(figFRCY,'papersize',ps,'paperposition',[0 0 ps])
set(figADYmean,'papersize',ps,'paperposition',[0 0 ps])
set(figADYeddy,'papersize',ps,'paperposition',[0 0 ps])
%
figHs  = figure('units','centimeters');
figHs.Position(3:4)  = ps;
set(figHs ,'papersize',ps,'paperposition',[0 0 ps])
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
% get mean and eddy velocities for advection terms
fileInfo = ncinfo(rotFile);
variableNames = {fileInfo.Variables.Name};
if ~ismember('Urot_mean',variableNames)
    U = mean(ncread(momFile,'umean'),3);
    V = mean(ncread(momFile,'vmean'),3);
    [~,U,V,~,~,~]=get_vel_decomposition_reGRID(Umean,Vmean,info.dx,info.dy);
else
    U = ncread(rotFile,'Urot_mean');
    V = ncread(rotFile,'Vrot_mean');
end
u     = ncread(rotFile,'Urot');
v     = ncread(rotFile,'Vrot');
% $$$ eta = ncread(rotFile,'eta');
% $$$ eta = eta-ETAmean;
%
disp('using time-averaged waterlevel... bug in source code')
ETA     = mean(ncread(momFile,'etamean'),3);
%
%% create depth mask (min-depth-resolved=0.01m, min-depth-normalize=0.1m)
h         = dep+ETA;% dep+eta;
mask      = h>0.01;
h(~mask)  = 0;
H         = dep+ETA;
MASK      = H>0.01;
H         = max(H, 0.1);
%
%% Wave Momentum Flux Terms:
UU = U.*U.*MASK;
VV = V.*V.*MASK;
UV = U.*V.*MASK;
uu = mean( u.*u.*h, 3, 'omitnan' ).*MASK./H;
uv = mean( u.*v.*h, 3, 'omitnan' ).*MASK./H;
vv = mean( V.*V.*h, 3, 'omitnan' ).*MASK./H;
%
% small scale filter
Nflty = 10/info.dy; if ~mod(Nflty,2), Nflty=Nflty+1;, end
Nfltx = 10/info.dx; if ~mod(Nfltx,2), Nfltx=Nfltx+1;, end
flt  = hanning(Nflty)*hanning(Nfltx)'; flt = flt./sum(flt(:));
uu = conv2(uu,flt,'same');
uv = conv2(uv,flt,'same');
vv = conv2(vv,flt,'same');
%
UU = conv2(UU,flt,'same');
UV = conv2(UV,flt,'same');
VV = conv2(VV,flt,'same');
%
% calculate eddy momentum terms
tmp = 0*uu;
tmp(:,2:end-1) = (uu(:,[3:end]) - uu(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uu(:,[2,end]) - uu(:,[1,end-1]))./(1*info.dx);

duudx = tmp;

tmp = 0*uv;
tmp(2:end-1,:) = (uv([3:end],:) - uv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (uv([2,end],:) - uv([1,end-1],:))./(1*info.dy);

duvdy = tmp;

tmp = 0*vv;
tmp(2:end-1,:) = (vv([3:end],:) - vv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (vv([2,end],:) - vv([1,end-1],:))./(1*info.dy);

dvvdy = tmp;

tmp = 0*uv;
tmp(:,2:end-1) = (uv(:,[3:end]) - uv(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uv(:,[2,end]) - uv(:,[1,end-1]))./(1*info.dx);

duvdx = tmp;

% Calculate mean momentum terms:
tmp = 0*UU;
tmp(:,2:end-1) = (UU(:,[3:end]) - UU(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UU(:,[2,end]) - UU(:,[1,end-1]))./(1*info.dx);

dUUdx = tmp;

tmp = 0*UV;
tmp(2:end-1,:) = (UV([3:end],:) - UV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (UV([2,end],:) - UV([1,end-1],:))./(1*info.dy);

dUVdy = tmp;

tmp = 0*VV;
tmp(2:end-1,:) = (VV([3:end],:) - VV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (VV([2,end],:) - VV([1,end-1],:))./(1*info.dy);

dVVdy = tmp;

tmp = 0*UV;
tmp(:,2:end-1) = (UV(:,[3:end]) - UV(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UV(:,[2,end]) - UV(:,[1,end-1]))./(1*info.dx);

dUVdx = tmp;
%
ADXeddy = duudx + duvdy;
ADXmean = dUUdx + dUVdy;
ADYeddy = dvvdy + duvdx;
ADYmean = dVVdy + dUVdx;

%
%% Pressure Gradient and Mean Advection:
g = 9.8;
ETA  = conv2(ETA,flt,'same');
[PgrdY,PgrdX] = gradientDG(g*ETA);
PgrdY = -PgrdY./info.dy;
PgrdX = -PgrdX./info.dx;
%
%% Wave Forcing:
DxSxx = ncread(momFile,'DxSxx');
DySxy = ncread(momFile,'DySxy');
DySyy = ncread(momFile,'DySyy');
DxSxy = ncread(momFile,'DxSxy');
%
BrkDissX = ncread(momFile,'BrkDissX');
BrkDissY = ncread(momFile,'BrkDissY');
%
RSX = -mean(DxSxx,3,'omitnan').*MASK./H - mean(DySxy,3,'omitnan').*MASK./H + mean(BrkDissX,3,'omitnan').*MASK./H;
RSY = -mean(DySyy,3,'omitnan').*MASK./H - mean(DxSxy,3,'omitnan').*MASK./H + mean(BrkDissY,3,'omitnan').*MASK./H;
%
%% Friction:
FRCX = ncread(momFile,'FRCX');
FRCY = ncread(momFile,'FRCY');
%
FRCX = -mean(FRCX,3,'omitnan').*MASK./H;
FRCY = -mean(FRCY,3,'omitnan').*MASK./H;
%
%
%% Wave Height:
Hs = ncread(momFile,'Hsig');
Hs = mean(Hs,3);
if spread==0
    fltWave = hamming(25); fltWave = fltWave'/sum(fltWave);
    Hs = conv2(Hs,fltWave,'same');
end
%
xlims = [0 450]/(info.xc-50);
ylims = 0.5*ph/pw*450/info.lc*[-1 1];
%
%% Pressure gradient:
figure(figPGX)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,PgrdX), caxis(1e-2*[-1 1]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%% Wave Forcing:
figure(figRSX)
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,RSX), caxis(1e-2*[-1 1]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
pos = get(gca,'position');
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');    
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%% Friction:
figure(figFRCX)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,FRCX), caxis(1e-3*[-1 1]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%
%% Mean Advection
figure(figADXmean)
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,ADXmean), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');    
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%
%% Eddy Advection:
figure(figADXeddy)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,ADXeddy), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
[~,h] = contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5);
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%% Pressure gradient:
figure(figPGY)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,PgrdY), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%% Wave Forcing:
figure(figRSY)
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,RSY), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
pos = get(gca,'position');
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');    
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%% Friction:
figure(figFRCY)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,FRCY), caxis(1e-3*[-1 1]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%
%% Mean Advection
figure(figADYmean)
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,ADYmean), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');    
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%
%% Eddy Advection:
figure(figADYeddy)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,ADYeddy), caxis(1e-2*[-0.5 0.5]), colormap(cmocean('balance'))
hold on,
[~,h] = contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5);
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%
%
%% Wave Height
figure(figHs)
eval(sprintf('pos = ppos%d;',ii))
axes('units','centimeters','position',pos);
imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,Hs), caxis([0.5 1.2]*height(ii)), colormap(cmocean('thermal'))
hold on,
[~,h] = contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5);
set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
str = split(cbttl,'[');
title([str{1},'$=~',cblbl{ii},'$~[',str{2}],'fontsize',8)
if ii==1
    ylabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
else
    set(gca,'yticklabel',[])
    if ii==round(N/2)
        xlabel('$(x-x_{xl})/X_c$ [~]','interpreter','latex')
    end
end
%


end



figure(figPGX);
cm    = cmocean('balance');
clims = [-1 1];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'X-Pres. Grad.',' [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')

figname = [figDIR,'PressureGradX_',NAME,'.png'];
exportgraphics(figPGX,figname)

figure(figRSX);
cm    = cmocean('balance');
clims = [-1 1];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'X-Wave Frc.'; '[m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'WaveForceX_',NAME,'.png'];
exportgraphics(figRSX,figname)
%%
figure(figFRCX);
cm    = cmocean('balance');
clims = [-1 1];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'X-Friction'; '  [m/s$^2$]$\times 10^{-3}$'},'interpreter','latex')
figname = [figDIR,'FrictionX_',NAME,'.png'];
exportgraphics(figFRCX,figname)
%
figure(figADXmean);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'X-Mean Adv.';'  [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'AdvectionX_mean_',NAME,'.png'];
exportgraphics(figADXmean,figname)
%
%
figure(figADXeddy);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'X-Eddy Adv.';'  [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'AdvectionX_eddy_',NAME,'.png'];
exportgraphics(figADXeddy,figname)


figure(figPGY);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'Y-Pres. Grad.','[m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')

figname = [figDIR,'PressureGradY_',NAME,'.png'];
exportgraphics(figPGY,figname)

figure(figRSY);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'Y-Wave Frc.'; '  [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'WaveForceY_',NAME,'.png'];
exportgraphics(figRSY,figname)
%%
figure(figFRCY);
cm    = cmocean('balance');
clims = 1e-3*[-1 1];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'Y-Friction'; '  [m/s$^2$]$\times 10^{-3}$'},'interpreter','latex')
figname = [figDIR,'FrictionY_',NAME,'.png'];
exportgraphics(figFRCY,figname)
%
figure(figADYmean);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'Y-Mean Adv.';'  [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'AdvectionY_mean_',NAME,'.png'];
exportgraphics(figADYmean,figname)
%
%
figure(figADYeddy);
cm    = cmocean('balance');
clims = [-0.5 0.5];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
xlabel({'Y-Eddy Adv.';'  [m/s$^2$]$\times 10^{-2}$'},'interpreter','latex')
figname = [figDIR,'AdvectionY_eddy_',NAME,'.png'];
exportgraphics(figADYeddy,figname)
%
%
%
figure(figHs);
cm    = cmocean('thermal');
clims = [0.5 1.2];
cvals = clims(1):diff(clims)/255:clims(2);
cb = axes('units','centimeters','position',cbpos);
imagesc(0,cvals,reshape(cm,256,1,3))
set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
       'fontsize',6,'tickdir','out')
cb.Position(3)=0.75*cb.Position(3);
xlabel({'$H_\mathrm{s}/H_0$ [~]'},'interpreter','latex')
figname = [figDIR,'Hs_',NAME,'.pdf'];
exportgraphics(figHs,figname)
%
close all
end