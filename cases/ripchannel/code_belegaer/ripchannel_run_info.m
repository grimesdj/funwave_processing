function info = ripchannel_run_info(run,runWAVES,init,grid);
% 
% run  = {'planar','barred','terraced','LYYWXX'};
% runWAVES  = {'h05s05','h05s10','h05s20','h10s05','h10s10',h10s20,'h15s05','h15s10','h15s20'};
% on belegaer see: funwave_run_info.m or darwin_local_run_info.m for other test simulations
% 
if ~exist('init','var')
    init=0;
end

if ~exist('grid','var')
    grid = run;
% $$$ else
% $$$     cmsROOT = ['/scratch/grimesdj/ripchannel/',run,filesep];
% $$$     lclROOT = ['/data2/ripchannel/',run,filesep];
end
cmsROOT = '/scratch/grimesdj/ripchannel/';
lclROOT = '/data2/ripchannel/';


infoFile = sprintf([lclROOT,run,filesep,'mat_data/ripchannel_run_info_%s_%s.mat'],grid,runWAVES);

if exist(infoFile,'file') & ~init
    info = load(infoFile);
else
info.rootMOD = [cmsROOT,run,filesep];
info.rootOut = [info.rootMOD,filesep,grid,'_',runWAVES,filesep,'output',filesep];
info.rootDAT = [lclROOT,run,filesep];
info.rootSim = [info.rootDAT,filesep,grid,'_',runWAVES,filesep];
info.rootMat = [info.rootDAT,filesep,'mat_data',filesep];
info.rootInp = [info.rootDAT,'inputs',filesep];
info.rootName= ['funwave_',grid,'_',runWAVES,'_'];
info.fileName = infoFile;
info.runName  = [grid,'_',runWAVES];
switch grid
  case 'planar2D'
    info.slope = 0.03;
    info.Ly    = 3e3;
    info.dx    = 1.0;
    info.dy    = 1.0;
  case 'planar1D'
    info.slope = 0.03;
    info.Ly    = 3.0;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'barred1D'
    info.slope = 0.03;
    info.xc    = 212.0;
    info.wc    = 27.5;
    info.ac    = 2.26;
    info.Ly    = 3.0;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'terraced1D'
    info.slope = 0.03;
    info.xc    = 205.0;
    info.wc    = 36.6;
    info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'terraced1Ddx150'
    info.slope = 0.03;
    info.xc    = 205.0;
    info.wc    = 36.6;
    info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 1.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'terraced1Ddx050'
    info.slope = 0.03;
    info.xc    = 205.0;
    info.wc    = 36.6;
    info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'terraced1Ddx025'
    info.slope = 0.03;
    info.xc    = 205.0;
    info.wc    = 36.6;
    info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 0.25;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'planar1Ddx050'
    info.slope = 0.03;
% $$$     info.xc    = 205.0;
% $$$     info.wc    = 36.6;
% $$$     info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'planar1Ddx025'
    info.slope = 0.03;
% $$$     info.xc    = 205.0;
% $$$     info.wc    = 36.6;
% $$$     info.ac    = 1.81;
    info.Ly    = 3.0;
    info.dx    = 0.25;
    info.dy    = 1.0;
    info.is1D  = true;
  case 'plnr2D'
    info.slope = 0.03;
    info.Ly    = 3e3;    
    info.dx    = 0.5;
    info.dy    = 1.0;
  case 'bar2D'
    info.slope = 0.03;
    info.xc    = 212.0;
    info.wc    = 27.5;
    info.ac    = 2.26;
    info.Ly    = 3e3;
    info.dx    = 0.5;
    info.dy    = 1.0;
  case 'ter2D'
    info.slope = 0.03;
    info.xc    = 205.0;
    info.wc    = 36.6;
    info.ac    = 1.81;
    info.Ly    = 3e3;
    info.dx    = 0.5;
    info.dy    = 1.0;
  case {'barRip0','barRip1','barRip2','barRip3','bar1Rip0','bar1Rip1','bar1Rip2','bar1Rip3'}
    info.slope = 0.03;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.Ly    = 3e3;
    info.rc    = 1;
    if str2num(grid(4))==1
        info.xc    = 150.0;
        info.wc    = 16.83;
        info.ac    = 1.39;
    else
        info.xc    = 212.0;
        info.wc    = 27.5;
        info.ac    = 2.26;
    end
    %
    if str2num(grid(end))==0
        info.lc    = 100;
    elseif str2num(grid(end))==1
        info.lc    = 50;
    elseif str2num(grid(end))==2
        info.lc    = 25;
    elseif str2num(grid(end))==3
        info.lc    = 150;
    end
  case {'terRip0','terRip1','terRip2','terRip3','ter1Rip0','ter1Rip1','ter1Rip2','ter1Rip3'}
    info.slope = 0.03;
    info.dx    = 0.5;
    info.dy    = 1.0;
    info.Ly    = 3e3;
    info.rc    = 1;
    if str2num(grid(4))==1
        info.ac    = 1.11;
        info.xc    = 145.0;
        info.wc    = round(100*info.ac/info.slope*exp(-1/2))/100;
    else
        info.xc    = 205.0;
        info.ac    = 1.81;
        info.wc    = 36.6;
    end
    %
    if str2num(grid(end))==0
        info.lc    = 100;
    elseif str2num(grid(end))==1
        info.lc    = 50;
    elseif str2num(grid(end))==2
        info.lc    = 25;
    elseif str2num(grid(end))==3
        info.lc    = 150;
    end
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

% $$$ tmp = info;
% $$$ rootMOD = tmp.rootMOD;
% $$$ idx = strfind(rootMOD,':');
% $$$ if ~isempty(idx)
% $$$     rootMOD = rootMOD(idx+1:end);
% $$$ end
% $$$ tmp.rootMOD  = rootMOD;
% $$$ tmp.rootSim  = [rootMOD,filesep,run,'_',runWAVES,filesep];
% $$$ tmp.rootOut  = [tmp.rootSim,filesep,'output',filesep];
% $$$ tmp.rootDAT  = [rootMOD,filesep];
% $$$ tmp.rootMat  = [rootMOD,filesep,'mat_data',filesep];
% $$$ tmp.rootInp  = [rootMOD,filesep,'inputs',filesep];
% $$$ str          = sprintf('ripchannel_run_info_%s_%s.mat',run,runWAVES);
% $$$ tmp.fileName = [rootMOD,filesep,'mat_data',filesep,str];
% $$$ save([info.rootMat,filesep,str],'-struct','tmp')


if ~exist(info.rootSim,'dir'),
    eval(['!mkdir -p ', info.rootSim,'/output/']),
    
end

if ~exist(info.rootInp,'dir'),
    eval(['!mkdir -p ', info.rootInp]),
    
end

save(info.fileName,'-struct','info')