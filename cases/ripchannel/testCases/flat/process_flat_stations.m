clear all
close all

addpath(genpath('~/git/funwave/code/'))
%% compare wave statistics between 3 runs with different wavemaker discretizations:
% run0 = standard -60:60
% run1 = limit to +/-6 STDEV (over sampled)
% run2 = use resolution at peak frequency... limit to 6-STDEV
for run = 0:2
    %% First create the "info" structure that has directory and file information for code:
    info = struct([]);
    info(1).runName   = sprintf('run%d',run);
    info.rootName  = info.runName;
    info.gaugeFile = '/scratch/grimesdj/ripchannel/flat/gauge.txt';
    info.bathyFile = '/scratch/grimesdj/ripchannel/flat/flatBathy.mat';
    info.rootSim   = sprintf('/scratch/grimesdj/ripchannel/flat/output%d',run);
    info.rootMat   = '/scratch/grimesdj/ripchannel/flat/mat_data/';
    info.fileName  = [info.rootMat,filesep,info.rootName,'_info.mat'];

    %% Read input spectrum from "LOG.txt"
    fin = [info.rootSim,filesep,'LOG.txt'];
    [freq, dire, amp] = parse_wavemaker_info_from_LOG(fin,0);
    info.freq=freq;
    info.dire=dire;
    info.amp =amp;

    %% Read final Hs field (Mglob and Nglob from input file)
    Mglob = 1364;
    Nglob = 2999;
    HsFiles = dir([info.rootSim,filesep,'output',filesep,'Hsig_*']);
    info.Hs = 0;
    Nf      = length(HsFiles);
    for jj=1:Nf
        fin = [HsFiles(jj).folder,filesep,HsFiles(jj).name];
        fid = fopen(fin);
        dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
        info.Hs = info.Hs+dum';% so I'm transposing so rows are y-coord and columns are x-coord
        fclose(fid);
    end
    info.Hs = info.Hs/Nf;

    info = process_FUNWAVE_virtual_moorings(info);
    
    out(run+1)=info;
end

save('/scratch/grimesdj/ripchannel/flat/mat_data/flat_runs.mat')
