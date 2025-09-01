function [info,t0,dT0] = ripchannel_run_info_cmshpc(matDIR,runNAME)


    infoFile = [matDIR,filesep,'ripchannel_run_info_',runNAME,'.mat'];
    info     = load(infoFile);
    % update root directories to match rootMOD, etc.,
    info.rootDAT = info.rootMOD;
    info.rootSim = [info.rootDAT,info.runName,filesep];
    info.rootMat = [info.rootDAT,'mat_data',filesep];
    info.rootInp = [info.rootDAT,'inputs',filesep];
    % update paths in filenames
    fileName = split(info.fileName,'/');
    fileName(cellfun(@isempty,fileName))=[];
    info.fileName= [info.rootMat,fileName{end}];
    grid = split(info.runName,'_');
    info.bathyFile = [info.rootMat,grid{1},'_depth.mat'];
    info.gaugeFile = [info.rootMat,grid{1},'_gauge.txt'];
    
    % 3) load the output times and dts
    info.timeFile = [info.rootSim,'time_dt.out'];
    Tdt = load(info.timeFile);
    t0   = Tdt(:,1);
    dt0  = gradient(t0);
    dT0  = Tdt(:,2); clear Tdt
    info.dt = mean(dt0);
    save(info.fileName,'-struct','info')

    
end
