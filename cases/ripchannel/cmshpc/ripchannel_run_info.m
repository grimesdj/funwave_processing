function info = ripchannel_run_info(runBATHY,runWAVES);
% 
% runBATHY  = {'planar','barred','terraced','LYYWXX'};
% runWAVES  = {'h05s05','h05s10','h05s20','h10s05','h10s10',h10s20,'h15s05','h15s10','h15s20'};
% on belegaer see: funwave_run_info.m or darwin_local_run_info.m for other test simulations
% 
infoFile = sprintf('/data2/ripchannel/mat_data/ripchannel_run_info_%s_%s.mat',runBATHY,runWAVES);
if exist(infoFile,'file')
    info = load(infoFile);
else
info.rootMOD = ['cmshpc:/scratch/grimesdj/ripchannel/',runBATHY];
info.rootOut = [info.rootMOD,filesep,runBATHY,'_',runWAVES,filesep,'output',filesep];
info.rootDAT = ['/data2/ripchannel/',runBATHY,filesep];
info.rootSim = [info.rootDAT,filesep,runBATHY,'_',runWAVES,filesep];
info.rootMat = [info.rootDAT,filesep,'mat_data',filesep];
info.rootInp = [info.rootDAT,'inputs',filesep];
info.rootName= ['funwave_',runBATHY,'_',runWAVES,'_'];
info.fileName = infoFile;
info.runName  = [runBATHY,'_',runWAVES];
switch runBATHY
  case 'planar2D'
    info.slope = 0.03;
    info.Ly    = 3e3;
    info.dx    = 1.0;
    info.dy    = 1.0;
  case 'planar1D'
    info.slope = 0.03;
    info.Ly    = 3.0;
    info.dx    = 1.0;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'barred1D'
    info.slope = 0.03;
    info.Ly    = 3.0;
    info.dx    = 1.0;
    info.dy    = 1.0;
    info.is1D  = true;
end    
% parse the wave height/perios/dir/spread info
Hs = regexp(runWAVES,'(?<=h)(..)','match');
info.Hs = str2num(Hs{1})/10;
Tp = regexp(runWAVES,'(?<=t)(..)','match');
info.Tp = str2num(Tp{1});        
Dp = regexp(runWAVES,'(?<=d)(..)','match');
info.Dp = str2num(Dp{1});        
spread = regexp(runWAVES,'(?<=s)(..)','match');
info.spread=str2num(spread{1});
end

if ~exist(info.rootMat,'dir'),
    eval(['!mkdir -p ', info.rootMat]),
end

tmp = info;
rootMOD = tmp.rootMOD;
idx = strfind(rootMOD,':');
rootMOD = rootMOD(idx+1:end);
tmp.rootDAT  = [rootMOD,filesep,runBATHY,filesep];
tmp.rootSim  = [rootMOD,filesep,runBATHY,'_',runWAVES,filesep];
tmp.rootMat  = [rootMOD,filesep,'mat_data',filesep];
tmp.rootInp  = [rootMOD,filesep,'inputs',filesep];
str          = sprintf('ripchannel_run_info_%s_%s.mat',runBATHY,runWAVES);
tmp.fileName = [rootMOD,filesep,'mat_data',filesep,str];
save([info.rootMat,filesep,str],'-struct','tmp')


if ~exist(info.rootSim,'dir'),
    eval(['!mkdir -p ', info.rootSim,'/output/']),
    
end

if ~exist(info.rootInp,'dir'),
    eval(['!mkdir -p ', info.rootInp]),
    
end

save(info.fileName,'-struct','info')