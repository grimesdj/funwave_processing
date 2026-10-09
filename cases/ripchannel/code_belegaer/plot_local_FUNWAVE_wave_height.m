clear all
% close all
addpath('~/git/funwave/code/')
%
runBATHY = 'spreadRip';
% smooth wave height distribution
Nflt = 21; flt = hanning(Nflt); flt = flt./sum(flt);
%
Hall = [];
for spread = [0, 2, 4, 10, 20]
runID    = sprintf('barRip0_h10t10s%02dd00',spread);
%
info = prep_local_ripchannel_info(runBATHY,runID);
%
%
waveStatsFile = dir([info.rootMat,'*',runID,'*wave_statistics.nc']);
info.waveStatsFile = [waveStatsFile(1).folder,filesep,waveStatsFile(1).name];
%
x = ncread(info.waveStatsFile,'x');
y = ncread(info.waveStatsFile,'y');
%
depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
depFile = [depFile(1).folder,filesep,depFile(1).name];
h       = ncread(depFile,'dep');
%
ETA    = ncread(info.waveStatsFile,'eta');
H      = h+ETA;
mask   = (h+ETA)>0.1;
%
Hs = ncread(info.waveStatsFile,'Hs');
Hs = mean(Hs,1,'omitnan')';
Hs = conv(Hs,flt,'same');
Hall = cat(2,Hall,Hs);
end
N = size(Hall,2);
cm = cmocean('thermal',N+1);
cm = cm(1:N,:);
%
figure,
colororder(cm);
plot(x,Hall)
%
% $$$ PSI    = ncread(info.waveStatsFile,'PSI');
% $$$ PSI(~mask) = nan;
% $$$ PSIavg = mean(PSI,3);
% $$$ PSI    = PSI-PSIavg;
% $$$ 
% $$$ cm = cmocean('balance');
% $$$ 
% $$$ figure, imagesc(y,x,PSIavg'), colormap(cm), colorbar, caxis([-0.0118    0.0118]),set(gca,'ydir','normal')
