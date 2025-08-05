% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYlist = {'planar2D'};
reproc  = 1;% 1=reprocess ascii to mat
rmfiles = 0;% 1=remove original ascii files when finished
for ii=1:length(runBATHYlist)
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
for jj = 1:Ndirs
    %
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
    % need local virtual mooring location file
    gaugeFile = split(info.gaugeFile,'/');
    info.stationFile = [info.rootInp,gaugeFile{end}];
    %
    if reproc & jj>1
        % 3) load the output times and dts
        info.timeFile = [info.rootSim,'time_dt.out'];
        Tdt = load(info.timeFile);
        t   = Tdt(:,1);
        dt  = gradient(t);
        dT  = Tdt(:,2); clear Tdt
        info.dt = mean(dt);
        info.bathyFile = [info.rootMat,filesep,runBATHY,'_depth.mat'];
        save(info.fileName,'-struct','info')
        %
        % 4) convert the funwave output ascii files to .mat
        fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
    end
    %
    info.subDomain = [1 3000 1 500];
    %
    % 5) process full run wave stats
    fprintf('\t estimating wave stats 1hr simulation: %s \n', runID)    
    info = estimate_FUNWAVE_run_wave_stats(info);
    %
    % 6) estimate velocity decomposition (U_rot, U_irr), alongshore spectra, ...
    fprintf('\t estimating velocity deocmposition: %s \n', runID)        
    info = estimate_FUNWAVE_velocity_helmholtz_decomposition(info);
    info = process_FUNWAVE_virtual_moorings(info);    
    info = estimate_FUNWAVE_rotational_velocity_alongshore_spectra(info);
    info = estimate_FUNWAVE_spectral_energy_flux(info);    
    
    % 7) estimate rotational power input and rotational impulse
    fprintf('\t estimating power and impulse: %s \n', runID)            
    info = estimate_FUNWAVE_rotational_power_input(info);
    info = estimate_FUNWAVE_rotational_impulse(info);


end


end
