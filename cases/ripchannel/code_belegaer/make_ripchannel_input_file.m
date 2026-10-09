function make_ripchannel_input_file(info)
% uses info structure to set the grid and wave specs
% and example stable input file: /data2/ripchannel/funwave_input.sample
% to create an "input.txt" file for a funwave-tvd run.
%
ini  = [info.rootSim,filesep,'input.txt'];
fid0 = fopen(ini,'w');
fprintf(fid0,'TITLE=%s\n',info.runName);
%
BATHY   = regexp(info.runName,'^.+?(?=_)','match');
BATHY   = BATHY{1};
rootSIM = regexp(info.rootMOD,'(?<=:).*$','match');
if isempty(rootSIM);
    rootSIM = info.rootMOD;
else
    rootSIM = rootSIM{1};
end
fprintf(fid0,'DEPTH_FILE=%s\n',['../inputs',filesep,BATHY,'_depth.txt']);
fprintf(fid0,'Mglob=%d\n',info.Nx-1);
fprintf(fid0,'Nglob=%d\n',info.Ny-1);
fprintf(fid0,'DX=%1.2f\n',info.dx);
fprintf(fid0,'DY=%1.2f\n',info.dy);
fprintf(fid0,'Xc_WK=%4.1f\n',info.xWM);
fprintf(fid0,'DEP_WK=%2.1f\n',info.hWM);
if isfield('info','is1D')
    scale = sqrt(2);
else
    scale = 1;
end
fprintf(fid0,'Hmo=%2.2f\n',info.Hs/scale);
fprintf(fid0,'FreqPeak=%2.4f\n',1/info.Tp);
fprintf(fid0,'ThetaPeak=%2.4f\n',info.Dp);
fprintf(fid0,'Sigma_Theta=%2.4f\n',info.spread);
fprintf(fid0,'NumberStations=%d\n',info.Ng);
fprintf(fid0,'STATIONS_FILE=%s\n',['../inputs',filesep,BATHY,'_gauge.txt']);
fclose(fid0)
if isfield(info,'is1D')
    eval(['!cat /data2/ripchannel/funwave_input1D.sample >> ',ini])
elseif info.spread==0
    eval(['!cat /data2/ripchannel/funwave_input_zero_spread.sample >> ',ini])
else
    eval(['!cat /data2/ripchannel/funwave_input.sample >> ',ini])
end
