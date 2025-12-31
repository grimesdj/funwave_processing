%% grimes edited to skip reprocessing the first run... %%
 
% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'testRip'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 0;% 1=remove original ascii files when finished
recalc  = 1;% 1=recalculate run statistics
plot_all_momentum=0;
%
%
%
% $$$ for ii=1:length(runBATHYlist)
ii=1
runBATHY = runBATHYlist{ii};
%
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
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
% before beginning, need to specify which variable we're processing:
% for all cases archive but don't delete the following vars:
% $$$ vars = {'dep','eta','u','v','mask','BrkSrcX','BrkSrcY'};
%
% loop over run_dirs and process WaveAveragedFields
Ndirs  = length(run_dirs);
for jj = 1:Ndirs
% 1) get current run subdirectory to process:
runID    = run_dirs{jj};
fprintf('\n processing wave-averaged: %s - %s \n', runBATHY,runID)    
% 2) get the archived info structure:
infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
if length(infoFile)>1
    fprintf('\tmultiple run-info files for:\t %s\n',runID)
    fprintf('\tusing filename:\t\t\t %s\n',infoFile(1).name);
end
info  = load([infoFile(1).folder,filesep,infoFile(1).name]);
%
%
if reproc & jj>2
    info = prep_info_structure(info);
    %
    % run specific output grid info:
    spanx = 1;
    spany = 1;
    rngx  = [1 info.subDomain(4)+1];% rngx  = [1 info.Lx/info.dx];
    rngy  = [1 info.subDomain(2)+1];% rngy  = [1 info.Ly/info.dy];
    %
    % construct time vector for wave-averaged variables
    vars = {'dep','etawavg','uwavg','vwavg'};
% $$$     waveAvgFiles = dir([info.rootOut,filesep,'etawavg*']);
% $$$     NwaveAvg = length(waveAvgFiles);
    NwaveAvg = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_wavg);
    dt_lp = info.T_INTV_wavg*ones(NwaveAvg,1);
    t_lp = [1:NwaveAvg]*dt_lp(1);
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
    %
    % construct time vector for Radiation Stress variables
    vars = {'dep','BrkDissX','BrkDissY','DxSxx','DxSxy','DxUUH','DxUVH','DySxy','DySyy','DyUVH','DyVVH','FRCX','FRCY','PgrdX','PgrdY','Sxx','Syy','Sxy','umean','vmean','etamean'};
% $$$     BrkDissFiles = dir([info.rootOut,filesep,'BrkDissX*']);
% $$$     Nbrk = length(BrkDissFiles);
    Nbrk = floor((info.TOTAL_TIME-info.STEADY_TIME)/info.T_INTV_mean);
    dt_lp = info.T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    fLog = convert_funwave_output_to_single_NetCDF(info.rootOut,[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
end
%
%
if recalc & jj~=1
    info = estimate_FUNWAVE_run_statistics_WaveAvgVelocity(info);
end
%
%
% plot the run statistics
fout_stats    = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info)
fout_momentum = plot_FUNWAVE_run_momentum(info)
end
%
if plot_all_momentum
%% compare each case from the same "grid"...
grids = split(run_dirs,'_');
grids = grids(:,1);
%
%% 0) get a unique list of grids:
ugrid = unique(grids);
%
%
% 1) define figure parameters, create two figures (along, cross), and colormap
xm = 2;
ym = 2;
pw = 4;
ph = 2;
ag = 0.2;
%
ppos31 = [xm           ym           pw ph];
ppos32 = [xm+pw+ag     ym           pw ph];
ppos33 = [xm+2*(pw+ag) ym           pw ph];
ppos21 = [xm           ym+ph+ag     pw ph];
ppos22 = [xm+pw+ag     ym+ph+ag     pw ph];
ppos23 = [xm+2*(pw+ag) ym+ph+ag     pw ph];
ppos11 = [xm           ym+2*(ph+ag) pw ph];
ppos12 = [xm+pw+ag     ym+2*(ph+ag) pw ph];
ppos13 = [xm+2*(pw+ag) ym+2*(ph+ag) pw ph];
cbpos  = [xm+3*(pw+ag)+xm ym 2*ag 1.5*ph];
ps     = [2*xm+3*pw+8*ag 2*ym+3*ph+3*ag];
%
fig = figure('units','centimeters');
pos = get(fig,'Position');
set(fig,'Position',[pos(1:2) ps], 'Papersize',ps,'PaperPosition',[0 0 ps])
for ii=1:3,
    for jj=1:3
        eval(['ax',num2str(ii),num2str(jj),'=axes(''units'',''centimeters'',''position'',ppos',num2str(ii),num2str(jj),');'])
    end
end
%
% 2) loop over ugrid
for ii = 1:length(ugrid);
    % 3) find all runs with current ugrid
    idx = ismember(grids,ugrid(ii));
    % 4) loop over all runs on this grid
    for jj=1:length(idx)
        % 5) load the alongshore transects of x- and y-momentum
        infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
        info  = load([infoFile(1).folder,filesep,infoFile(1).name]);
        momFile = [info.rootMat,info.rootName,'MomentumTerms'];
        %
        % 4.2.1) estimate time-averages of:
        %        cross-shore terms:
        %             advection, 
        tmp1 = ncread(momFile,'DxUUH');
        tmp2 = ncread(momFile,'DyUVH');
        ADX  = mean(tmp1,3,'omitnan')+mean(tmp2,3,'omitnan'); clear tmp1 tmp2
        %             pressure grad,
        PGX  = ncread(momFile,'PgrdX');
        PGX  = mean(PgrdX,3,'omitnan');
        %             radiation stress+BrkDissX,
        tmp1 = ncread(momFile,'DxSxx');        
        tmp2 = ncread(momFile,'DySxy');
        tmp3 = ncread(momFile,'BrkDissX');        
        RSX  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan') - mean(tmp3,3,'omitnan');
        clear tmp1 tmp2 tmp3
        %        along-shore terms:
        %             advection, 
        tmp1 = ncread(momFile,'DyVVH');
        tmp2 = ncread(momFile,'DxUVH');
        ADY  = mean(tmp1,3,'omitnan'); clear tmp1 tmp2
        %             pressure grad,
        PGY  = ncread(momFile,'PgrdY');
        PGY  = mean(PGY,3,'omitnan');
        %             radiation,
        tmp1 = ncread(momFile,'DySyy');
        tmp2 = ncread(momFile,'DxSxy');
        tmp3 = ncread(momFile,'BrkDissY');
        RSY  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan') - mean(tmp3,3,'omitnan');
        clear tmp1 tmp2 tmp3
        %
        % 6) add curve to figure, add value to list for colorbar
        
    end
end
% 7) place colorbar, axis labels, etc.
%
%
end
%
%
%% fast time wave files
Ndirs  = length(run_dirs);
for jj = 1:Ndirs
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
    % 3) load the output times and dts
    info.timeFile = [info.rootSim,'time_dt.out'];
    Tdt = load(info.timeFile);
    t0   = Tdt(:,1);
    dt0  = gradient(t0);
    dT0  = Tdt(:,2); clear Tdt
    info.dt = mean(dt0);
% $$$     save(info.fileName,'-struct','info')
    %
    % 4) convert the funwave output ascii files to .mat
    vars = {'dep','eta','mask','BrkSrcX','BrkSrcY'};
    fLog = convert_funwave_output_to_NetCDF(info.rootOut,[info.rootMat,info.rootName],vars,t0,dT0,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
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
