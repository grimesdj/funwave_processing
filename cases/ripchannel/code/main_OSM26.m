% some workflow is done on belegaer:/data2/ripchannel/code/

% first, run the main_process_####.m for all runs = #####,
% then copy the "_dep.nc", "_MomentumTerms.nc", and "_velocity_decomposition.nc" files to belegaer.

% then add the details for run #### to the pre-processing portion of:
compile_WaveAvgVelocity_stats


% then modify and run:
compare_WaveAvgVelocity_stats

