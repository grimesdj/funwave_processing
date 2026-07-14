clear all
close all
%%
% Make runfiles for 1D flat runs. Trying to get rid of standing waves in domain.
inputFile = '/storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/flat1D';
rootDir   = '/scratch/grimesdj/ripchannel/flat1D/'

% Need to decide on a parameter space to use:
%
% 1) begin with 50-m width, increase width by ~2.5% for each run,
% 2) decrease the magnitude of R_sponge and CDsponge with each iteration...

% So far, the best result was: L=51.25, Rrate = 0.8718, CDsponge = 5.

L0 = 70; Lrate = 1.00;
R_sponge = 0.8718; Rrate = 1.00;
CDsponge = 7;      Crate = 1.025;


for iter = 0:15

    % build info structure:
    info.runName = ['run', num2str(iter)];
    info.rootSim = [rootDir, filesep, 'run', num2str(iter)];
    info.rootMat = [rootDir, filesep, 'mat_data'];
    info.fileName= [info.rootMat, filesep, info.runName, '_info.mat'];
    info.rootOut = [info.rootSim,filesep,'output',filesep];
    info.SpongeWidth = Lrate.^iter*L0;
    info.Rsponge     = Rrate.^(-iter)*R_sponge;
    info.CDsponge    = Crate.^(-iter)*CDsponge;    

    if ~exist(info.rootSim,'dir')
        eval(['!mkdir -p ',info.rootSim])
    end
    
    % make input file:
    ini  = [info.rootSim,filesep,'input.txt'];
    fid0 = fopen(ini,'w');
    fprintf(fid0,'TITLE=%s\n',info.runName);
    fprintf(fid0,'Sponge_west_width=%f\n',info.SpongeWidth);
    fprintf(fid0,'Sponge_east_width=%f\n',info.SpongeWidth);        
    fprintf(fid0,'R_sponge=%f\n',info.Rsponge);
    fprintf(fid0,'CDsponge=%f\n',info.CDsponge);    
    fclose(fid0);

    % add the input fields that are not changing   
    eval(['!cat /storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/flat1D/input.template >> ',ini])    

    % add the run_name to the slurm file-list
    slurmRunList{iter+1} = info.runName;
end

sub   = [rootDir,filesep,'flat1D.slurm'];
fid1  = fopen(sub,'w');
fprintf(fid1,'#!/bin/bash\n');
fprintf(fid1,'# Array of run directories to submit\n');
fprintf(fid1,'run_dirs=(');
fprintf(fid1,'"%s" ',slurmRunList{:});
fprintf(fid1,')\n');

eval(['!cat /storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/flat1D/slurm_batch1D.sample >> ', sub]);
eval(['!chmod 777 ',sub])

