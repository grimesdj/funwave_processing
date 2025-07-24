% code to be launched on cms-hpc "cuttlefish"
% 1) requires the input bathymetry name as top-dir
% 1.1) need to loop over this if necessary...
runBATHY = 'planar2D';
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
reproc  = 1;
rmfiles = 0;
vars = {'dep','eta','u','v','mask','nubrk','p','q'};
%
% loop over run_dirs
Ndirs  = length(run_dirs);
for jj = 1:Ndirs
    %
    % 1) get current run subdirectory to process:
    runID    = run_dirs{jj};
    % 2) get the archived info structure:
    infoFile = dir([runBathy,'*','info','*',runID,'.mat']);
    if length(infoFile)>1
        fprintf('multiple run-info files for: %s\n',runID)
        fprintf('using filename: %s\n',infoFile(1).name);
    end
    info  = load([infoFile(1).folder,filesep,infoFile(1).name]);
    %
    % 3) load the output times and dts
    info.timeFile = [info.rootOut,'time_dt.out'];
    save(info.fileName,'-struct','info')
    Tdt = load(info.timeFile);
    t   = Tdt(:,1);
    dt  = gradient(t);
    dT  = Tdt(:,2); clear Tdt
    %
    % 4) convert the funwave output ascii files to .mat
    fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)



end
