% currently working on runs in:
runDay = 'OSM2026';
% for all cases archive but don't delete the following vars:
rmfiles = 0;
vars = {'dep','eta','u','v','mask','nubrk','p','q'};
%
runIDs = {'sig10uni','sig20','sig20uni'}; % {'full1mres','sig05','sig05uni','sig10'};
for ii=1:length(runIDs)
runID = runIDs{ii};
info  = darwin_run_info(runDay,runID);
info
%
% load the output times and dts
info.timeFile = [info.rootOut,'time_dt.out'];
save(info.fileName,'-struct','info')
Tdt = load(info.timeFile);
t   = Tdt(:,1);
dt  = gradient(t);
dT  = Tdt(:,2); clear Tdt
%
fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
%
% 
first30sec_time = [info.first30sOut,filesep,'time_dt.out'];
Tdt = load(first30sec_time);
t   = Tdt(:,1);
dt  = gradient(t);
dT  = Tdt(:,2); clear Tdt
fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
end
%
%
% $$$ %%%%%%%% first process the full runs: %%%%%%%%%%%%%%%
% $$$ runID = 'full';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ runID = 'narrow';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ %
% $$$ runID = 'uniform';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ %
% $$$ %%%%%%%% second process the first 30s of each run: %%%%%%%%%%%%%%%
% $$$ runID = 'full';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ runID = 'uniform';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ runID = 'narrow';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ runID = 'uniform_narrow';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
%
% $$$ %
% $$$ runID = 'narrow_normal';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ runID = 'narrow_normal'
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ runID = 'uniform_narrow'
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
% $$$ %
% $$$ %
% $$$ %
% $$$ runID = 'uniform_narrow_normal'
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.first30sOut,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
% $$$ %
% $$$ fLog = convert_funwave_output_to_mat(info.first30sOut,[info.rootMat,info.rootName,'first30sec_'],vars,t,dT,1,1,rmfiles,300)
%
% $$$ runID = 'uniform_narrow_normal';
% $$$ info  = darwin_run_info(runDay,runID);
% $$$ %
% $$$ % load the output times and dts
% $$$ info.timeFile = [info.rootSim,'time_dt.out'];
% $$$ save(info.fileName,'info')
% $$$ Tdt = load(info.timeFile);
% $$$ t   = Tdt(:,1);
% $$$ dt  = gradient(t);
% $$$ dT  = Tdt(:,2); clear Tdt
%
% $$$ fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
%
