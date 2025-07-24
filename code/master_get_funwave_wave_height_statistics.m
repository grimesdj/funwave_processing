% $$$ %
% $$$ % get run-directory info
% $$$ root = '/Users/derekgrimes/funwave/sz2Dturb/';
% $$$ bath = 'plnr02';
% $$$ runs  = {'sig04','sig10','sig20','sig40'};
% $$$ Hs  = 0.5;
% $$$ Tp  = 8;
% $$$ %
% $$$ plotter=0;
% $$$ %
% $$$ for ii=1:length(runs)
% $$$ run = runs{ii};
% $$$ % get model Hs and surf-zone width
% $$$ [Hs_x,x,eta_bar,h_bar,sig_eta,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(root,bath,run,Hs,Tp,plotter);
% $$$ %
% $$$ fout = [root,'mat_data/WaveHeight_',bath,'_',run,'.mat'];
% $$$ save(fout,'Hs_x','x','eta_bar','h_bar','xsl','xsz','Lsz')
% $$$ %
% $$$ end

% $$$ % code to estimate fun--wave stats
if ~exist('runDay','var') | ~exist('runID','var')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0929 run00                  %%%
%%                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root   = '/home/derek/projects/ShortCrests/MOD/funwave/';
runDay = '0929';
runID  = 'run00';
disp('processing default case: 0929 run00')
end
[rootSim,rootOut,rootMat,rootName,bathyFile,t,dt,x,dx,y,dy,h,Tp,Hs] = funwave_run_info(runDay,runID);
%
disp(['working on: ', runDay, ' ', runID])
%
nx = length(x);
ny = length(y);
SubDomain = [1 ny 1 nx];
%
[Hs_x,x,eta_bar,h_bar,sig_eta,Snn_xy,Snn_x,Snn_wg,freqs,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(rootMat,rootName,bathyFile,Hs,Tp);
% $$$ [Hs_x,x,eta_bar,h_bar,sig_eta,freqs,Snn_xy,Snn_x,Snn_wg,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(rootMat,rootName,bathyFile,run,Hs,Tp);
%
