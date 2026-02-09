function info = prep_local_ripchannel_info(runBATHY,runID);
%
% USAGE: info = prep_local_ripchannel_info(runBATHY,runID);
%
% Example: runBATHY = 'spreadRip';
%          runID    = 'barRip0_h10t10s10d00';

rootDIR  = '/data2/ripchannel/';
runDIR   = [rootDIR,runBATHY];
matDIR   = [runDIR,filesep,'mat_data',filesep];
%
infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
%
%
info.rootMOD = runDIR;
info.rootMat = matDIR;
