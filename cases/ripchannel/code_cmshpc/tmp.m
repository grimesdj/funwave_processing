%% 2) loop over the simulations,
iter = 0;
for cc = 1:length(CASES)
    CASE = CASES{cc};
    runDIR   = [rootDIR,CASE];
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
    Ndirs = length(run_dirs);
    %% 3) loop over each run in simulation,






        for rr = 1:length(run_dirs);
        iter = iter+1;
        RUN  = run_dirs{rr};
        %% 4) load info file:
        fprintf('\n processing: %s %s \n', CASE,RUN)    
        infoFile = dir([matDIR,filesep,'*','info','*',RUN,'.mat']);
        info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
        INFO(iter)=info;
        
        %% 4.1) get grid & wave parameters for run:
        grid_waves = split(info.runName,'_');
        grid  = grid_waves{1};
        waves = split(grid_waves{2},{'h','t','s','d'});
        %
        height(iter)    = str2num(waves{2})/10;
        period(iter)    = str2num(waves{3});
        spread(iter)    = str2num(waves{4});
        direction(iter) = str2num(waves{5});
        %
        bar_location(iter)            = info.xc;
        bar_width(iter)               = info.wc;
        bar_amplitude(iter)           = info.ac;
        channel_length(iter)          = info.lc;
        channel_amplitude_ratio(iter) = info.rc;

        %% 4.1.1) define label sizes/shapes/colors?
       
        %% 4.2) load post-processed variables:
        rotFile = [info.rootMat,'funwave_',RUN,'_velocity_decomposition.nc'];
        % get (x,y,dep)
        x   = ncread(rotFile,'x');
        y   = ncread(rotFile,'y');
        dep = ncread(rotFile,'dep');
        domain_size(iter,:) = size(dep);
        fprintf('\n Domain Size: Nx = %d Ny = %d \n',domain_size(iter,2),domain_size(iter,1))
        end
end    
