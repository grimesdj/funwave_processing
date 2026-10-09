function make_ripchannel_slurm_file(info,run_dirs)
%
% create the slurm job script
% BATHY   = regexp(info.runName,'^.+?(?=_)','match');

[RUN,mtch] = strsplit(info.rootDAT,filesep);
RUN(cellfun(@isempty,RUN))=[];
RUN   = RUN{end};
sub  = [info.rootDAT,filesep,RUN,'.slurm'];
fid1 = fopen(sub,'w');
fprintf(fid1,'#!/bin/bash\n');
fprintf(fid1,'# Array of run directories to submit\n');
fprintf(fid1,'run_dirs=(');
fprintf(fid1,'"%s" ',run_dirs{:});
fprintf(fid1,')\n');
if isfield(info,'is1D')
    eval(['!cat /data2/ripchannel/slurm_batch1D.sample >> ', sub]);
else    
    eval(['!cat /data2/ripchannel/slurm_wait.sample >> ', sub]);
end
eval(['!chmod 777 ',sub])