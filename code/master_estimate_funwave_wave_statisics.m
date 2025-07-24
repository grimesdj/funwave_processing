%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0929 run00                  %%%
%%                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root   = '/home/derek/projects/ShortCrests/MOD/funwave/';
runDay = '0929';
runID  = 'run00';
%
[rootSim,rootOut,rootMat,rootName,bathyFile,t,dt,x,dx,y,dy,h] = funwave_run_info(runDay,runID);
%
disp(['working on: ', runDay, ' ', runID])
%
[Hs_x,x,eta_bar,h_bar,sig_eta,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics([rootMat,rootName],bathyFile,run,Hs,Tp,subDomain,plotter);
