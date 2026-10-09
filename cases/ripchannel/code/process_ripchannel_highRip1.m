%% This is designed to process runs developed during summer 2026 (after OSM26)
%
% The new set of cases are: {spreadRip, gapRip, gapRip1, highRip, highRip1}
%
% The simulations are in: /scratch/grimesdj/ripchannel/
% The primary differences are:
% simulation duration = 6500s after 2500s spin-up (~30 time larger than repeat period of neighboring modes 1/df)
% sponge-layer: (direct) R_sponge=0.8718, (friction) Cdsponge = 2.8844,
% momentum term integration time: 6500 s

% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
rootDIR = '/scratch/grimesdj/ripchannel/';

% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'highRip1'};
reproc  = 1;% 1=reprocess ascii/binary to netcdf
rmfiles = 1;% 1=remove original ascii/binary files when finished
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
% pre-OSM26 issue: disp('ONLY PROCESSING NON-ZERO SPREAD CASES [2:5,7:10,12:15]')
for jj = 1:Ndirs
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
    % construct time vector for wave-averaged variables
    vars = {'dep','etawavg','uwavg','vwavg'};
    NwaveAvg = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_wavg);
    dt_lp = info.T_INTV_wavg*ones(NwaveAvg,1);
    t_lp = [1:NwaveAvg]*dt_lp(1);
    % binary output files... faster archiving
    isBINARY = 1;
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t_lp,dt_lp,info.dx,info.spanx,info.rngx,info.dy,info.spany,info.rngy,rmfiles,300,[-inf inf],isBINARY,info.Nx-1,info.Ny-1);
    % construct time vector for Radiation Stress variables
    vars = {'BrkDissX','BrkDissY','DxSxx','DxSxy','DxUUH','DxUVH','DySxy','DySyy','DyUVH','DyVVH','FRCX','FRCY','PgrdX','PgrdY','Sxx','Syy','Sxy','umean','vmean','etamean','Hsig'};
    Nbrk = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_mean);
    dt_lp = info.T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    fLog = convert_funwave_output_to_single_NetCDF(info.rootOut,[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,info.spanx,info.rngx,info.dy,info.spany,info.rngy,rmfiles,300,[-inf inf],isBINARY,info.Nx-1,info.Ny-1);
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
close all
end
end
end

