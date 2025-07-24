function [rootSim,rootOut,rootMat,rootName,bathyFile,t,dt,x,dx,y,dy,h,Tp,Hs,Dir,Sig,Gam,xWM,hWM,angle,dateTime] = funwave_run_info(runDay,runID);

rootMOD = '/home/derek/projects/ShortCrests/MOD/funwave/';
rootDAT = '/data0/ShortCrests/MOD/';
rootSim = [rootMOD,runDay,filesep,runID,filesep];
rootMat = [rootDAT,'mat_data',filesep];
rootOut = [rootDAT,runDay,filesep,runID,filesep,'output',filesep];
rootName= ['funwave_',runDay,'_',runID,'_'];

switch runDay
  case '0929'
    switch runID
      case 'run00'
        bathyFile = [rootMat,'bathy_small_from_20130905.mat'];
        Tp = 1/0.1;
        Hs = 1.51;
        Dir=0;
        Sig=0;
        Gam=1.9;
        xWM = [];% I don't remember
        hWM = [];
        % input file used gam=1, Hs=1.85 and fp=0.09 to model Hs=1.58, fp=0.1, 
        % but model Hs was a bit too big (and had too much LF energy).
        %
        % trying again w/ gam=1.4, Hs=1.77, and fp=0.1, much better
        % still too much LF stuff shallow, spectum much broader too...
        %
        % v2: trying gam=1.9 and fp=0.10125
        % trying gam=1.9 and fp=0.10125
        %
        % v3: trying gam=2.2, fp = 0.1
        %
        % v4: trying gam=1.9, fp = 0.1, froude = 5
      case 'run1D'
        bathyFile = [rootMat,'bathy_1D_50sLine_from_20130905.mat'];
        Tp = 1/0.0825;
        Hs = 1.3;
        Dir=0;
        Sig=0;
        Gam=3.3;
        xWM=1123.0;
        hWM=9;
        angle = 2;% grid rotation angle from FRF-coords
        dateTime='09291200';
        % similar to run00, but improved bathy 
        % trying:  gam=1.9, Hmo=1.78, fp = 0.1, froude = 3
        % retrying:
        %--First run: 07/03/2023
        %   Waterlevel=0;
        %   Hs=1.78      -->significant wave height too big
        %   GammaTMA=1.9 -->spectrum too broad
        %--Second run: 07/05/2023
        %   Waterlevel=0.36; based on FRF tide-gauge
        %   Hs= 1.51     --> based on AWAC Hs; too much energy at HF
        %   GammaTMA=3.3 --> use JoNSWAP     ; surf-zone spectrum too broad
        %   FreqPeak=0.1 --> based on AWAC Tp; in SZ mod Tp=0.11, obs Tp=0.093
        %--Third run: 07/06/2023
        %   Hs= 1.5      --> 
        %   GammaTMA=5.0 --> 
        %   FreqPeak=0.093->
        %--Fourth run: 07/07 --> still too much energy at swell and lower frequencies!
        %   Hs=1.45
        %   GammaTMA=3.3
        %   FreqPeak=0.09
        %--Fifth run: 08/07 --> peak was smaller than observed and at 9.2 seconds
        %   Hs=1.4
        %   GammaTMA=1.1;
        %   FreqPeak=0.09; 
        %--Sixth run: 12/22/23
        %   Hs=1.4
        %   GammaTMA=3.3;
        %   FreqPeak=0.085;
        %   cut-off freq: 1/20 1/5
        % --seventh run: 12/27
        %   changed:
        %   Hs=1.35
        %   FreqPeak=0.0625;
        % --eigth run: 12/28
        %   GammaTMA=3.3;
        %   FreqPeak=0.085;
        %   cutoff: 1/20 and 1/5
        %   water-level=0.18
        %   changed:
        %     Hs=1.3
        %     FreqMin=0.0625;
        % --ninth run: 12/28
        %   GammaTMA=1.9;
      case 'run1Dv2_hWM9m'
        bathyFile = [rootMat,'bathy_1D_50sLine_from_20130905_hWM9m.mat'];
        Tp = 10;
        Hs = 1.45;
        Dir=0;
        Sig=0;
        Gam=0;
        xWM=1123.0;
        hWM=9;
        angle = 2;% grid rotation angle from FRF-coords
        dateTime='09291200';
        % identical to run1D, except the wave amplitude spectrum is provided from the shoaled AWAC 1D spectrum.
      case 'run1D_awac'
        bathyFile = [rootMat,'bathy_1D_50sLine_from_20130905_hWM7m.mat'];
        Tp = 10;
        Hs = 1.00;
        Dir=0;
        Sig=0;
        Gam=0;
        xWM=692;
        hWM=7;
        angle = 2;% grid rotation angle from FRF-coords
        dateTime='09291200';
        % identical to run1D, except the wave amplitude spectrum is provided from the shoaled AWAC 1D spectrum.
      case 'run1D_wr'
        bathyFile = [rootMat,'bathy_1D_50sLine_from_20130905_hWM7m.mat'];
        Tp = 10;
        Hs = 1.00;
        Dir=0;
        Sig=0;
        Gam=0;
        xWM=692;
        hWM=7;
        angle = 2;% grid rotation angle from FRF-coords
        dateTime='09291200';
        % identical to run1D, except the wave amplitude spectrum is provided from the shoaled WR 1D spectrum.
      case 'run000'
        bathyFile = [rootMat,'bathy_small_from_20131023.mat'];
        Tp = 1/0.10125;
        Hs = 1.51;
        Dir=0;
        Sig=0;
        Gam=1.9;
      case 'run01'
        % full domain simulation based on the FA2023 rotated/periodized bathymetry
        % follows after run1D where I was tinkering to match Snn between obs/mod
        % along the 50's ADV instrument line, including the ROD.
        bathyFile = [rootMat,'bathy_large_from_20130905_hWM9m.mat'];
        %        rootOut   = ['/data0/ShortCrests/MOD/0929/run01/output_old_v2/'];
        Tp = 10;
        Hs = 1.3;
        Dir=0;
        Sig=20;
        Gam=1.9;
        xWM=1123.0;
        hWM=9.0;
        angle = 2;% grid rotation angle from FRF-coords
        dateTime='09291200';
      case 'run02'
        % this is a repeat of run01 with shallower WM and using AWAC data
        % it blows up due to non-linearity a/h~0.2 at wavemaker
        bathyFile = [rootMat,'bathy_large_from_20130905_hWM7m.mat'];
        %        
        Tp = 10;
        Hs = 1.5;
        Dir=0;
        Sig=0;
        Gam=0;
        xWM=692;
        hWM=7;
        angle=2;
        dateTime='09291200';
      case 'run03'
        % this is a repeat of run01 with shallower WM and using WR data
        % same as run02, blows up!
        bathyFile = [rootMat,'bathy_large_from_20130905_hWM7m.mat'];
        %        
        Tp = 10;
        Hs = 1.5;
        Dir=0;
        Sig=0;
        Gam=0;
        xWM=692;
        hWM=7;
        angle=2;
        dateTime='09291200';
    end
end

timeFile = [rootSim,'/','time_dt.out'];
if ~exist(timeFile,'file')
    timeFile = [rootOut,'time_dt.out'];
end

Tdt = load(timeFile);
if isempty(Tdt)
    t=[];
    dt=[];
    dT=[];
else
    t = Tdt(:,1);
    dt= gradient(t);
    dT= Tdt(:,2); clear Tdt
end

load(bathyFile,'x','y','h')
dx = x(2)-x(1);
dy = y(2)-y(1);

% wave maker is measured in meters from SW most cell
xWM = xWM+x(1);