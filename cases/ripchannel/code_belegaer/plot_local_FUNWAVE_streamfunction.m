clear all
% close all
addpath('~/git/funwave/code/')
%
runBATHY = 'spreadRip';
runID    = 'barRip0_h10t10s00d00';
%
info = prep_local_ripchannel_info(runBATHY,runID);
%
%
rotVelFile = dir([info.rootMat,'*',runID,'*velocity_decomposition.nc']);
info.rotVelFile = [rotVelFile(1).folder,filesep,rotVelFile(1).name];

x = ncread(info.rotVelFile,'x');
y = ncread(info.rotVelFile,'y');

depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
depFile = [depFile(1).folder,filesep,depFile(1).name];
h       = ncread(depFile,'dep');

ETA    = ncread(info.rotVelFile,'eta');
H      = h+ETA;
Havg   = mean(H,3,'omitnan');
mask   = (h+ETA)>0.1;

VORT    = ncread(info.rotVelFile,'VORT');
VORT(~mask) = nan;

momFile = dir([info.rootMat,'*',runID,'*MomentumTerms.nc']);
info.momFile = [momFile(1).folder,filesep,momFile(1).name];
Umean = ncread(info.momFile,'umean'); Umean = mean(Umean,3,'omitnan');
Vmean = ncread(info.momFile,'vmean'); Vmean = mean(Vmean,3,'omitnan');

[~,Vx] = gradientDG(Vmean/info.dx);
[Uy,~] = gradientDG(Umean/info.dy);
VORTavg=Vx-Uy;

VORT    = VORT-VORTavg;

PSI    = ncread(info.rotVelFile,'PSI');
PSI(~mask) = nan;

%
[PSIavg,Urot_avg,Vrot_avg,~,~,~]=get_vel_decomposition_reGRID(Umean,Vmean,info.dx,info.dy);        
% PSIavg = mean(PSI,3);
PSI    = PSI-PSIavg;




cm = cmocean('balance');

figure, imagesc(y,x,PSIavg'), colormap(cm), colorbar, 
hold(gca,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
caxis([-10 10]),set(gca,'ydir','normal')


[~,Vavg] = gradientDG(PSIavg./info.dx);
[Uavg,~] = gradientDG(-PSIavg./info.dy);

figure, imagesc(y,x,Vavg'), colormap(cm), colorbar, 
hold(gca,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
caxis([-0.25 0.25]),set(gca,'ydir','normal')


% zoom into ripchannel:
iX = find(x>info.xc-2*info.wc);
iY = find(y>info.Ly/2-3*info.lc & y<info.Ly/2+3*info.lc);

% transport:
tmp  = Uavg(iY,iX).*Havg(iY,iX);
msk  = tmp<=0;
tmp(msk)=nan;
T    = sum(tmp,1,'omitnan');

tmp = Havg(iY,iX);
tmp(msk)=nan;
Urip = T./sum(tmp,1,'omitnan');