info   = run_info_MOD(runDay,runID);
%
info.timeFile = [info.rootSim,'time_dt.out'];
save(info.fileName,'info')
Tdt = load(info.timeFile);
t = Tdt(:,1);
dt= gradient(t);
dT= Tdt(:,2); clear Tdt
rmfiles = 0;
vars = {'dep','eta','u','v','mask','nubrk','p','q'};
%
fLog = convert_funwave_output_to_mat(info.rootOut,[info.rootMat,info.rootName],vars,t,dT,1,1,rmfiles,300)
