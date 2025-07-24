function info = darwin_run_info(runDay,runID);
% on darwin:
% runDay = 'AGU2024'
% runID  = {'full','uniform','narrow'};

% on belegaer see: funwave_run_info.m for other test simulations
infoFile = sprintf('/lustre/xg-ees240064/users/3800/AGU2024/mat_data/darwin_run_info_%s.mat',runID);
if exist(infoFile,'file')
    info = load(infoFile);
else
info.rootMOD = '/lustre/xg-ees240064/users/3800/';
info.rootDAT = info.rootMOD;
info.rootSim = [info.rootMOD,filesep,runDay,filesep,runID,filesep];
info.rootMat = [info.rootDAT,filesep,runDay,filesep,'mat_data',filesep];
info.rootOut = [info.rootDAT,runDay,filesep,runID,filesep,'output',filesep];
info.rootName= ['funwave_',runDay,'_',runID,'_'];
info.fileName = infoFile;
info.runName  = [runDay,runID];
switch runDay
  case 'AGU2024'
    switch runID
      case 'full'
        info.rootOut     = [info.rootDAT,runDay,filesep,runID,filesep,'output_full_bathy_and_waves',filesep];
        info.first30sOut = [info.rootDAT,runDay,filesep,runID,filesep,'first30seconds_full_bathy_and_waves',filesep];                
% $$$       waveSource= 'swn26wr'
% $$$              hWM= 9
% $$$          deltaWM= 0.500000000000000
% $$$            angle= 2
% $$$         dateTime= '09291200'
% $$$       waterLevel= 0.180000000000000
% $$$        bathyFile= '/data0/ShortCrests/MOD/mat_data//funwave_0929_swn0_bathy.mat'
% $$$              xWM= 1098
% $$$         rootSwan= '/home/derek/projects/ShortCrests/MOD/funwave/../swan/0929/1200/'
% $$$        swanXlims= [-2 16218]
% $$$        swanYlims= [-5000 7000]
% $$$        swanDelta= 10
% $$$         waveFile= '/home/derek/projects/ShortCrests/MOD/funwave/0929/swn0//waves.swn26wr'
% $$$         waveNdir= 988
% $$$         waveNfrq= 988
% $$$               Tp= 11.750881316098708
% $$$               Hs= 1.365655519498581
% $$$          freqRNG= [0.055560000000000 0.228330000000000]
% $$$          direRNG= [-50.323760000000007 51.397940000000006]
% $$$               dx= 0.500000000000000
% $$$               dy= 1
% $$$               Nx= 2551
% $$$               Ny= 1501
% $$$               Ly= 1500
% $$$               Lx= 1275
% $$$               wl= 0.180000000000000
      case 'narrow'
        info.rootOut     = [info.rootDAT,runDay,filesep,runID,filesep,'output_full_bathy_and_narrow_waves',filesep];
        info.first30sOut = [info.rootDAT,runDay,filesep,runID,filesep,'first30seconds_full_bathy_and_narrow_waves',filesep];
      case 'uniform'
        info.rootOut     = [info.rootDAT,runDay,filesep,runID,filesep,'output_uniform_bathy_and_full_waves',filesep];
        info.first30sOut = [info.rootDAT,runDay,filesep,runID,filesep,'first30seconds_uniform_bathy_and_full_waves',filesep];      
    end
% $$$   case '0929'
% $$$     switch runID
% $$$       case 'run1D'
% $$$         info.is1D    = 1;
% $$$         info.bathyFile = [info.rootMat,'bathy_1D_50sLine_from_20130905.mat'];
% $$$         info.Tp = 1/0.0825;
% $$$         info.Hs = 1.3;
% $$$         info.Dir=0;
% $$$         info.Sig=0;
% $$$         info.Gam=3.3;
% $$$         info.xWM=1123.0;
% $$$         info.hWM=9;
% $$$         info.angle = 2;% grid rotation angle from FRF-coords
% $$$         info.dateTime='09291200';
% $$$       case 'run1Dv2_hWM9m'
% $$$         info.bathyFile = [info.rootMat,'bathy_1D_50sLine_from_20130905_hWM9m.mat'];
% $$$         info.Tp = 10;
% $$$         info.Hs = 1.45;
% $$$         info.Dir=0;
% $$$         info.Sig=0;
% $$$         info.Gam=0;
% $$$         info.xWM=1123.0;
% $$$         info.hWM=9;
% $$$         info.angle = 2;% grid rotation angle from FRF-coords
% $$$         info.dateTime='09291200';
% $$$         % identical to run1D, except the wave amplitude spectrum is provided from the shoaled AWAC 1D spectrum.
% $$$       case 'run1D_awac'
% $$$         info.bathyFile = [info.rootMat,'bathy_1D_50sLine_from_20130905_hWM7m.mat'];
% $$$         info.Tp = 10;
% $$$         info.Hs = 1.00;
% $$$         info.xWM=692;
% $$$         info.hWM=7;
% $$$         info.angle = 2;% grid rotation angle from FRF-coords
% $$$         info.dateTime='09291200';
% $$$         % identical to run1D, except the wave amplitude spectrum is provided from the shoaled AWAC 1D spectrum.
% $$$       case 'run1D_wr'
% $$$         info.bathyFile = [info.rootMat,'bathy_1D_50sLine_from_20130905_hWM7m.mat'];
% $$$         info.Tp = 10;
% $$$         info.Hs = 1.00;
% $$$         info.xWM=692;
% $$$         info.hWM=7;
% $$$         info.angle = 2;% grid rotation angle from FRF-coords
% $$$         info.dateTime='09291200';
% $$$         % identical to run1D, except the wave amplitude spectrum is provided from the shoaled WR 1D spectrum.
% $$$       case 'run01'
% $$$         % full domain simulation based on the FA2023 rotated/periodized bathymetry
% $$$         % follows after run1D where I was tinkering to match Snn between obs/mod
% $$$         % along the 50's ADV instrument line, including the ROD.
% $$$         info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM9m.mat'];
% $$$         %        rootOut   = ['/data0/ShortCrests/MOD/0929/run01/output_old_v2/'];
% $$$         info.Tp = 10;
% $$$         info.Hs = 1.3;
% $$$         info.Dir=0;
% $$$         info.Sig=20;
% $$$         info.Gam=1.9;
% $$$         info.xWM=1123.0;
% $$$         info.hWM=9.0;
% $$$         info.deltaWM=2;
% $$$         info.angle = 2;% grid rotation angle from FRF-coords
% $$$         info.dateTime='09291200';
% $$$       case 'run0awac'
% $$$         % this is a repeat of run01 with shallower WM and using AWAC data
% $$$         % waves too non-linear a/h=0.2        
% $$$         info.FRFbathyDate = '20130905';
% $$$         info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM7m.mat'];
% $$$         info.waveSource='awac';
% $$$         info.xWM=692;
% $$$         info.hWM=7;
% $$$         info.deltaWM=1;
% $$$         info.angle=2;
% $$$         info.dateTime='09291200';
% $$$         info.waterLevel = 0.18;
% $$$       case 'run0wr'
% $$$         % this is a repeat of run01 with shallower WM and using WR data
% $$$         % waves too non-linear a/h=0.2
% $$$         info.FRFbathyDate = '20130905';        
% $$$         info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM7m.mat'];
% $$$         info.waveSource='wr';
% $$$         info.xWM=692;
% $$$         info.hWM=7;
% $$$         info.deltaWM=1;
% $$$         info.angle=2;
% $$$         info.dateTime='09291200';
% $$$         info.waterLevel = 0.18;
% $$$       case 'run1awac'
% $$$         % this is a repeat of run01 with 0905 bathy and 9m WM
% $$$         % a/h = 0.16
% $$$         info.FRFbathyFile = '/home/derek/projects/ShortCrests/MOD/data/FRF_geomorphology_elevationTransects_survey_20130905.nc';
% $$$         info.SKIbathyFile = '/home/derek/projects/ShortCrests/MOD/data/RODDebbieTimeLagFix_20131002.xyz';
% $$$         info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM7m.mat'];
% $$$         info.waveSource='awac';
% $$$         info.hWM=9;
% $$$         info.deltaWM=0.5;
% $$$         info.angle=2;
% $$$         info.dateTime='09291200';
% $$$         info.waterLevel = 0.18;
% $$$       case 'run1wr'
% $$$         % this is a repeat of run01 with 0905 bathy and 9m WM
% $$$         % a/h = 0.16
% $$$         info.FRFbathyFile = '/home/derek/projects/ShortCrests/MOD/data/FRF_geomorphology_elevationTransects_survey_20130905.nc';
% $$$         info.SKIbathyFile = '/home/derek/projects/ShortCrests/MOD/data/RODDebbieTimeLagFix_20131002.xyz';
% $$$         info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM9m.mat'];
% $$$         info.waveSource='wr'
% $$$         info.hWM=9;
% $$$         info.deltaWM=0.5;
% $$$         info.angle=2;
% $$$         info.dateTime='09291200';
% $$$         info.waterLevel = 0.18;
% $$$       case 'swn0'
% $$$         % this is a repeat of run01 with 0905 bathy and 9m WM
% $$$         % a/h = 0.16
% $$$         info.FRFbathyFile = '/home/derek/projects/ShortCrests/MOD/data/FRF_geomorphology_elevationTransects_survey_20130905.nc';
% $$$         info.SKIbathyFile = '/home/derek/projects/ShortCrests/MOD/data/RODDebbieTimeLagFix_20131002.xyz';
% $$$         %        info.bathyFile = [info.rootMat,'bathy_large_from_20130905_hWM9m.mat'];
% $$$         info.waveSource='swn26wr'
% $$$         info.hWM=9;
% $$$         info.wl=0.18;
% $$$         info.deltaWM=0.5;
% $$$         info.angle=2;
% $$$         info.dateTime='09291200';
% $$$         info.waterLevel = 0.18;
% $$$         %
% $$$         rootMOD= '/home/derek/projects/ShortCrests/MOD/funwave/';
% $$$          rootDAT= '/data0/ShortCrests/MOD/';
% $$$          rootSim= '/home/derek/projects/ShortCrests/MOD/funwave/0929/swn0/';
% $$$          rootMat= '/data0/ShortCrests/MOD/mat_data/';
% $$$          rootOut= '/data0/ShortCrests/MOD/0929/swn0/output/';
% $$$         rootName= 'funwave_0929_swn0_';
% $$$         fileName= '/home/derek/projects/ShortCrests/MOD/mat_data/funwave_run_info_swn0.mat';
% $$$          runName= '0929swn0';
% $$$     FRFbathyFile= '/home/derek/projects/ShortCrests/MOD/data/FRF_geomorphology_elevationTransects_survey_20130905.nc';
% $$$     SKIbathyFile= '/home/derek/projects/ShortCrests/MOD/data/RODDebbieTimeLagFix_20131002.xyz';
% $$$       waveSource= 'swn26wr';
% $$$              hWM= 9;
% $$$          deltaWM= 0.5000;
% $$$            angle= 2;
% $$$         dateTime= '09291200';
% $$$       waterLevel= 0.1800;
% $$$        bathyFile= '/data0/ShortCrests/MOD/mat_data//funwave_0929_swn0_bathy.mat';
% $$$              xWM= 1098;
% $$$         rootSwan= '/home/derek/projects/ShortCrests/MOD/funwave/../swan/0929/1200/';
% $$$        swanXlims= [-2 16218];
% $$$        swanYlims= [-5000 7000];
% $$$        swanDelta= 10;
% $$$         waveFile= '/home/derek/projects/ShortCrests/MOD/funwave/0929/swn0//waves_with_boundIG.swn26wr';
% $$$         waveNdir= 901;
% $$$         waveNfrq= 901;
% $$$               Tp= 11.7509;
% $$$               Hs= 1.2885;
% $$$          freqRNG= [0.0556 0.2306];
% $$$          direRNG= [-49.9327 49.8033];
% $$$               dx= 0.5000;
% $$$               dy= 1;
% $$$               Nx= 2551;
% $$$               Ny= 1501;
% $$$               Ly= 1500;
% $$$               Lx= 1275;
% $$$     end
end
end

save(info.fileName,'-struct','info')


if ~exist(info.rootSim,'dir')
    disp('making the simulation directory')
    eval(['!mkdir -p ',info.rootSim])
end


if ~exist(info.rootOut,'dir')
    disp('making the output directory and symbolic link')
    eval(['!mkdir -p ',info.rootOut])
    eval(['!ln -s ', info.rootOut, ' ', [info.rootSim,filesep,'output' ]])
end
