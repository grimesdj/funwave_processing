clear all
close all
%
build = 0;
if build
    % first define a set of runs
    RUN = {'gapRip1'};%{'spreadRip'};%{'barred1D', 'terraced1D'};% 
    run_dirs = build_ripchannel_run_files(RUN);
end
%
%% Fall 2025: issues with Hrms vs Hsig in 1D and 2D %%
% $$$ % explore the issue with 1D wave heights being too large.
% $$$ get_model_1D_Hm0_correction_factor

%% NSF 2026 Report Bathy %%
run = 'spreadRip'
[~,~,~,~,runWAVESlist,grids] = ripchannel_parameter_space(run)
%
reinitialize = 0;
runWAVES     = runWAVESlist{1};
grid         = grids{1};
% info = ripchannel_run_info(run,runWAVES,reinitialize,grid)
info = prep_local_ripchannel_info(run,[grid,'_',runWAVES]);
info = make_ripchannel_bathy(info,1)
%%

return
%% OSM26 Figures CODE %%
%% build/plot bathymetry to show rip-channel structure
run = 'spreadRip'
[~,~,~,~,runWAVESlist,grids] = ripchannel_parameter_space(run)
%
reinitialize = 0;
runWAVES     = runWAVESlist{1};
grid         = grids{3};
% info = ripchannel_run_info(run,runWAVES,reinitialize,grid)
info = prep_local_ripchannel_info(run,[grid,'_',runWAVES]);

%% END OSM26 %%

