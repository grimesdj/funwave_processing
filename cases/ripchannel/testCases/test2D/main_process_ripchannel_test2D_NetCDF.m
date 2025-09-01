% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'test2D'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 0;% 1=remove original ascii files when finished
%
% need info from input file
T_INTV_mean = 200;
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
% loop over run_dirs
Ndirs  = length(run_dirs);
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
if reproc & jj>1
    % update root directories to match rootMOD, etc.,
    info.rootDAT = info.rootMOD;
    info.rootSim = [info.rootDAT,info.runName,filesep];
    info.rootMat = [info.rootDAT,'mat_data',filesep];
    info.rootInp = [info.rootDAT,'inputs',filesep];
    % update paths in filenames
    fileName = split(info.fileName,'/');
    fileName(cellfun(@isempty,fileName))=[];
    info.fileName= [info.rootMat,fileName{end}];
    grid = split(info.runName,'_');
    info.bathyFile = [info.rootMat,grid{1},'_depth.mat'];
    info.gaugeFile = [info.rootMat,grid{1},'_gauge.txt'];
    %
    % define subDomain for analysis
    info.subDomain = [1 info.Ny-1 1 round(400/info.dx)-1];
    %
    % 3) load the output times and dts
    info.timeFile = [info.rootSim,'time_dt.out'];
    Tdt = load(info.timeFile);
    t0   = Tdt(:,1);
    dt0  = gradient(t0);
    dT0  = Tdt(:,2); clear Tdt
    info.dt = mean(dt0);
    save(info.fileName,'-struct','info')
    %
    spanx = 1;
    spany = 1;
    %
    rngx  = [1 400/info.dx];
    rngy  = [1 info.Ly/info.dy];
    % 4) convert the funwave output ascii files to .mat
    rootOut = split(info.rootOut,':');
    vars = {'dep','eta','u','v','p','q','mask','BrkSrcX','BrkSrcY'};
    fLog = convert_funwave_output_to_NetCDF(rootOut{end},[info.rootMat,info.rootName],vars,t0,dT0,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
    %
    % construct time vector for Radiation Stress variables
    vars = {'BrkDissX','BrkDissY','DxSxx','DxSxy','DxUUH','DxUVH','DySxy','DySyy','DyUVH','DyVVH','FRCX','FRCY','PgrdX','PgrdY','Sxx','Syy','Sxy'};
    BrkDissFiles = dir([rootOut{end},filesep,'BrkDissX*']);
    Nbrk = length(BrkDissFiles);
    dt_lp = T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    fLog = convert_funwave_output_to_single_NetCDF(rootOut{end},[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,300);
end
%
%
info = estimate_FUNWAVE_run_statistics_BrkSrc(info);
end
% $$$ end
