%% error in input files (zero-spread are non-binary i/o) corrupted first run: barRip0_h05t10s00d00
%% processing s==0 here...


% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
rootDIR = '/scratch/grimesdj/ripchannel/';
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'highRip'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 1;% 1=remove original ascii files when finished
recalc  = 1;% 1=recalculate run statistics
plotter = 1;% 1=plot run statistics
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
disp('ONLY PROCESSING NON-ZERO SPREAD CASES 2:2:12')
for jj = [1:2:Ndirs]
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
    % binary output files... faster archiving
    isBINARY = 0;
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300,[-inf inf],isBINARY,info.Nx-1,info.Ny-1);
    % construct time vector for Radiation Stress variables
    vars = {'BrkDissX','BrkDissY','DxSxx','DxSxy','DxUUH','DxUVH','DySxy','DySyy','DyUVH','DyVVH','FRCX','FRCY','PgrdX','PgrdY','Sxx','Syy','Sxy','umean','vmean','etamean'};
    Nbrk = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_mean);
    dt_lp = info.T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    fLog = convert_funwave_output_to_single_NetCDF(info.rootOut,[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300,[-inf inf],isBINARY,info.Nx-1,info.Ny-1);
end
%
%
if recalc 
    info = estimate_FUNWAVE_run_statistics_WaveAvgVelocity(info);
end
%
%
% plot the run statistics
if plotter
fout_stats    = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info)
fout_momentum = plot_FUNWAVE_run_momentum(info)
fout_momentum_offline = plot_FUNWAVE_offline_momentum_budget(info)
end
end
end

