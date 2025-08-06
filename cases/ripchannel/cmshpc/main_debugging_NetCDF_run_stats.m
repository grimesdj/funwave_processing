% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'planar2D'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 0;% 1=remove original ascii files when finished
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
vars = {'dep','eta','u','v','mask','nubrk','p','q'};
%
% loop over run_dirs
Ndirs  = length(run_dirs);
% $$$ for jj = 1:Ndirs
jj=15
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
info.subDomain = [1 info.Ny-1 1 round(500/info.dx)-1];
%
if reproc 
    % 3) load the output times and dts
    info.timeFile = [info.rootSim,'time_dt.out'];
    Tdt = load(info.timeFile);
    t0   = Tdt(:,1);
    dt0  = gradient(t0);
    dT0  = Tdt(:,2); clear Tdt
    info.dt = mean(dt0);
    info.bathyFile = [info.rootMat,filesep,runBATHY,'_depth.mat'];
    save(info.fileName,'-struct','info')
    %
    spanx = 1;
    spany = 1;
    % 4) convert the funwave output ascii files to .mat
    rootOut = split(info.rootOut,':');
    fLog = convert_funwave_output_to_NetCDF(rootOut{end},[info.rootMat,info.rootName],vars,t0,dT0,info.dx,spanx,info.dy,spany,rmfiles,300)
end
%
%
info = estimate_FUNWAVE_run_statistics(info);
% $$$ end
% $$$ end
