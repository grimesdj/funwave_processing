function info = prep_info_structure(info);
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
    %
    % Use the input.txt file to populate T_INTV_MEAN, and T_INTV_WAVG
    fin = [info.rootSim, 'input.txt'];
    fid = fopen(fin);
    while ~feof(fid)
        line = fgetl(fid);
        str  = split(line,'=');
        if length(str)==1, continue, end
        var  = deblank(str{1});
        val  = str{2};
        switch var
          case {'TOTAL_TIME'}
            info.TOTAL_TIME = str2num(val);
          case {'STEADY_TIME'}
            info.STEADY_TIME = str2num(val);                        
          case {'T_INTV_mean'}
            info.T_INTV_mean = str2num(val);
          case {'T_INTV_wavg'}
            info.T_INTV_wavg = str2num(val);
        end
    end
    fclose(fid)
    % define subDomain for analysis... this is run specific!
    Lx = 500;
    fprintf('\nRestricting cross-shore analysis to x<=%f\n',Lx)
    info.Lx = Lx;
    info.subDomain = [1 info.Ny-1 1 round(Lx/info.dx)-1];
    info.spanx = 1;
    info.spany = 1;
    %
    % remove the remote hostname from output directory
    rootOut = split(info.rootOut,':');
    info.rootOut = rootOut{end};
    save(info.fileName,'-struct','info')    
end

