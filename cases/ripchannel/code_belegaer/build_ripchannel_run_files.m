function run_dirs = build_ripchannel_run_files(RUN);
%
% USAGE:
%
% builds a set of Funwave-TVD run folders and input grids for a pre-defined
% parameter space specified by:
% RUN:
%        'planar1D'
%        'planar2D'
%        'barred2D'
% loop over runs to generate all permutations
reinitialize=1;
for ii = 1:length(RUN);
    run_dirs = {};
    run      = RUN{ii};
    [HEIGHT,PERIOD,DIR,SPREAD,runWAVESlist,grids] = ripchannel_parameter_space(run);
    if isempty(grids)
        for jj=1:length(runWAVESlist)
            runWAVES = runWAVESlist{jj};
            %
            % create info structure
            info = ripchannel_run_info(run,runWAVES,reinitialize);
            %
            % generate bathymetry
            plotter = 0;
            if (jj) == 1
                plotter = 1;
            end
            info = make_ripchannel_bathy(info,plotter);
            %
            % construct funwave input file from info structure
            make_ripchannel_input_file(info);
            %
            % keep a running list of runs to zip
            % also need these to create slurm script
            run_dirs = cat(1,run_dirs,info.runName);
        end
    else
    for kk=1:length(grids)
        grid = grids{kk};
        for jj=1:length(runWAVESlist)
            runWAVES = runWAVESlist{jj};
            %
            % create info structure
            info = ripchannel_run_info(run,runWAVES,reinitialize,grid);
            %
            % generate bathymetry
            plotter = 0;
            if (jj) == 1
                plotter = 1;
            end
            info = make_ripchannel_bathy(info,plotter);
            %
            % construct funwave input file from info structure
            make_ripchannel_input_file(info);
            %
            % keep a running list of runs to zip
            % also need these to create slurm script
            run_dirs = cat(1,run_dirs,info.runName);
        end
    end
    end
        make_ripchannel_slurm_file(info,run_dirs)
        save([info.rootMat,'runs_to_process.mat'],'run_dirs')
        eval(['!tar --exclude="mat_data_TINTVmean200" -czvf ../',run,'.tar.gz   -C ',info.rootDAT,' .'])
end