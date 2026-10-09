clear all
close all

addpath(genpath('~/git/funwave/code/'))
%% compare wave statistics between 3 runs with different wavemaker discretizations:
runDirs = {'newSponge','planar1D'}
% $$$ run=0;
for run = 0:1
    %% First create the "info" structure that has directory and file information for code:
    info = struct([]);
    info(1).runName   = runDirs{run+1};
    info.rootName  = info.runName;
    info.rootSim   = ['/data2/ripchannel/planar1D/',info.runName];
    info.rootMat   =  '/data2/ripchannel/planar1D/newSponge/mat_data/';
    info.fileName  = [info.rootMat,filesep,info.rootName,'_info.mat'];

    if ~exist(info.rootMat,'dir')
        eval(['!mkdir -p ',info.rootMat])
    end
    
    %% Read input spectrum from "LOG.txt"
    fin = [info.rootSim,filesep,'LOG.txt'];
    [freq, dire, amp] = parse_wavemaker_info_from_LOG(fin,0);
    if isempty(freq)
        [freq, dire, amp] = parse_wavemaker_info_LOG_aggregated(fin,0);
        if isempty(freq)
            disp('no wavemaker informamtion read from LOG.txt')
        end
    end
                 
    info.freq=freq;
    info.dire=dire;
    info.amp =amp;

    %% Read final Hs field (Mglob and Nglob from input file)
    Mglob = 1364;
    Nglob = 2;
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

% $$$     info = process_FUNWAVE_virtual_moorings(info);
    save(info.fileName,'-struct','info')
    out(run+1)=info;
end

save('/data2/ripchannel/planar1D/newSponge/mat_data/planar_1D_runs.mat')
