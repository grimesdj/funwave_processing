%% Grimes edited to only process s00 runs --> lines 35 & 86: jj = [1, 6, 11]

% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
rootDIR = '/scratch/grimesdj/ripchannel/';
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'spreadRip'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 1;% 1=remove original ascii files when finished
recalc  = 1;% 1=recalculate run statistics 
%
%
for ii=1:length(runBATHYlist)
runBATHY = runBATHYlist{ii};
%
runDIR   = [rootDIR,runBATHY];
matDIR   = [runDIR,filesep,'mat_data'];
%
% the list of run directories are saved in:
load([matDIR,filesep,'runs_to_process.mat'])
% brings in cell array: run_dirs
% for example,
% run_dirs =
%   6x1 cell array
%    {'planar1D_h05t08s00d00'}
%    {'planar1D_h05t10s00d00'}
%    {'planar1D_h10t08s00d00'}
%    {'planar1D_h10t10s00d00'}
%    {'planar1D_h15t08s00d00'}
%    {'planar1D_h15t10s00d00'}
%
% loop over run_dirs
Ndirs  = length(run_dirs);
for jj = [1 6 11]%1:Ndirs
% 1) get current run subdirectory to process:
runID    = run_dirs{jj};
fprintf('\n processing: %s %s \n', runBATHY,runID)    
% 2) get the archived info structure:
infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
if length(infoFile)>1
    fprintf('\tmultiple run-info files for:\t %s\n',runID)
    fprintf('\tusing filename:\t\t\t %s\n',infoFile(1).name);
end
info  = load([infoFile(1).folder,filesep,infoFile(1).name]);
%
%
if reproc 
    info = prep_info_structure(info);
    %
    % run specific output grid info:
    spanx = 1;
    spany = 1;
    rngx  = [1 info.subDomain(4)+1];% rngx  = [1 info.Lx/info.dx];
    rngy  = [1 info.subDomain(2)+1];% rngy  = [1 info.Ly/info.dy];
    % construct time vector for wave-averaged variables
    vars = {'dep','etawavg','uwavg','vwavg'};
    NwaveAvg = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_wavg);
    dt_lp = info.T_INTV_wavg*ones(NwaveAvg,1);
    t_lp = [1:NwaveAvg]*dt_lp(1);
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
    % construct time vector for Radiation Stress variables
    vars = {'BrkDissX','BrkDissY','DxSxx','DxSxy','DxUUH','DxUVH','DySxy','DySyy','DyUVH','DyVVH','FRCX','FRCY','PgrdX','PgrdY','Sxx','Syy','Sxy','umean','vmean','etamean'};
    Nbrk = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_mean);
    dt_lp = info.T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    fLog = convert_funwave_output_to_single_NetCDF(info.rootOut,[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
end
%
%
if recalc 
    info = estimate_FUNWAVE_run_statistics_WaveAvgVelocity(info);
end
%
%
% plot the run statistics
fout_stats    = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info)
fout_momentum = plot_FUNWAVE_run_momentum(info)
end
end

% exit

% check for instantaneous files and proceed if they exist...
%% fast time wave files
for jj = [1 6 11]%1:Ndirs
% 1) get current run subdirectory to process:
runID    = run_dirs{jj};
fprintf('\n processing fast-time: %s - %s \n', runBATHY,runID)    
% 2) get the archived info structure:
infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
if length(infoFile)>1
    fprintf('\tmultiple run-info files for:\t %s\n',runID)
    fprintf('\tusing filename:\t\t\t %s\n',infoFile(1).name);
end
info  = load([infoFile(1).folder,filesep,infoFile(1).name]);
%
%
if reproc
    spanx = 1;
    spany = 1;
    rngx  = [1 info.subDomain(4)+1];% rngx  = [1 info.Lx/info.dx];
    rngy  = [1 info.subDomain(2)+1];% rngy  = [1 info.Ly/info.dy];
    % 3) load the output times and dts
    info.timeFile = [info.rootSim,'time_dt.out'];
    if ~exist(info.timeFile,'file')
        fprintf('No fast-time output for: %s-%s',runBATHY,runID)
        continue
    end
    Tdt = load(info.timeFile);
    t0   = Tdt(:,1);
    dt0  = gradient(t0);
    dT0  = Tdt(:,2); clear Tdt
    info.dt = mean(dt0);
% $$$     save(info.fileName,'-struct','info')
    %
    % 4) convert the funwave output ascii files to .mat
    vars = {'eta','mask','BrkSrcX','BrkSrcY'};
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t0,dT0,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
    save(info.fileName,'-struct','info')
end
%
if recalc 
    info = estimate_FUNWAVE_run_statistics_Waves(info);
end
%
%
% plot the run statistics
fout_stats    = plot_FUNWAVE_run_Wave_statistics(info)
end
