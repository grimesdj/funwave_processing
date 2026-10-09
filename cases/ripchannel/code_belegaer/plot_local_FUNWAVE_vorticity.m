clear all
% close all
addpath('~/git/funwave/code/')
%
runBATHY = 'spreadRipLong';
runID    = 'barRip0_h10t10s10d00';
%
info = prep_local_ripchannel_info(runBATHY,runID);
info.rotVelFile = [info.rootMat, info.rootName, 'velocity_decomposition.nc'];
%
info = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info);
%
% $$$ rotVelFile = dir([info.rootMat,'*',runID,'*velocity_decomposition.nc']);
% $$$ info.rotVelFile = [rotVelFile(1).folder,filesep,rotVelFile(1).name];
% $$$ 
% $$$ x = ncread(info.rotVelFile,'x');
% $$$ y = ncread(info.rotVelFile,'y');
% $$$ 
% $$$ depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
% $$$ depFile = [depFile(1).folder,filesep,depFile(1).name];
% $$$ h       = ncread(depFile,'dep');
% $$$ 
% $$$ ETA    = ncread(info.rotVelFile,'eta');
% $$$ H      = h+ETA;
% $$$ mask   = (h+ETA)>0.1;
% $$$ 
% $$$ VORT    = ncread(info.rotVelFile,'VORT');
% $$$ VORT(~mask) = nan;
% $$$ VORTavg = mean(VORT,3);
% $$$ VORT    = VORT-VORTavg;
% $$$ 
% $$$ 
% $$$ PSI    = ncread(info.rotVelFile,'PSI');
% $$$ PSI(~mask) = nan;
% $$$ PSIavg = mean(PSI,3);
% $$$ PSI    = PSI-PSIavg;
% $$$ 
% $$$ cm = cmocean('balance');
% $$$ 
% $$$ figure, imagesc(y,x,PSIavg'), colormap(cm), colorbar, caxis([-0.0118    0.0118]),set(gca,'ydir','normal')