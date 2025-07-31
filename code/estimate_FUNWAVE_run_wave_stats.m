function info = estimate_FUNWAVE_run_wave_stats(info);
%
% usage:  info = estimate_FUNWAVE_run_wave_stats(info);
%
% calculate wave frequency spectra, wave height, front statistics
% and breaking-wave rotational forcing. Analysis is limited to the region:
%       25 m <= x <=  400 m
%     -200 m <= y <= 1300 m
% get the input (x,y,h) grid file
load(info.bathyFile)
x0=x;
y0=y;
h0=h;
%
% full alongshore domain
if ~isfield(info,'subDomain')
    iX = find(x0>=0  & x0<=400);
    nx = length(iX);
    ny = length(y0);
    subDomain = [1 ny iX(1) iX(end)];
else
    subDomain = info.subDomain;
end
%
if ~all(ismember('first30sec',info.runName))
%
[Hs_xy,eta_bar,mask0,x,y,h,freq,Snn_xy,xsl] = calculate_funwave_wave_height_statistics_v2(info.rootMat,info.rootName,info.bathyFile,info.Hs,info.Tp,subDomain);
%
info.x_shoreline = xsl;
info.mask        = mask0;
info.eta_bar     = eta_bar;
%
fig0 = figure;
Hs_x = nanmean(Hs_xy,1);
p0 = plot(x,Hs_xy,'.k',x,Hs_x,'-r','markersize',1,'linewidth',1.5);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$H_s$ [m]','interpreter','latex')
set(gca,'xlim',[75 600],'ylim',[0 1.1*nanmax(Hs_xy(:))],'ticklabelinterpreter','latex','tickdir','out')
f0l1 = legend([p0(1),p0(end)]','$H_\mathrm{s}(x,y)$','$\bar{H}_\mathrm{s}(x)$');
set(f0l1,'location','southeast','interpreter','latex')
title(info.runName)
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ',[info.rootSim,filesep,'figures']])
end
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'Hs.pdf'])
%
end
%
% get front statistics
[mean_stats,binned_stats] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ [mL,sL,N,mLx,sLx,Lx_log_mean,Lx_log_std,Nx,xylog,Xbins] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% front_file = [info.rootMat,filesep,info.rootName,'bore_front_statistics.mat'];
%
% make a stats plot? what stats?
fig1 = figure;
p1 = plot(binned_stats.Xbins,binned_stats.Length,'-k',binned_stats.Xbins,exp( binned_stats.log_mean_length' + binned_stats.log_std_length'*[-1 1] ), '--r');
xlabel('crosshore [m]','interpreter','latex')
ylabel('$L$ [m]','interpreter','latex')
set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
f1l1 = legend([p1(1), p1(2)]','$\bar{L}(l)=\exp\left(\overline{\log(l)}\right)$','$\bar{L}\pm\mathrm{std}(L)$');
set(f1l1,'location','northeast','interpreter','latex')
title(info.runName)
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ',[info.rootSim,filesep,'figures']])
end
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'crest_length_vs_x.pdf'])
%
fig2 = figure;
p2 = plot(binned_stats.Xbins,binned_stats.N,'-k','linewidth',2);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$N$ [crests/frame]','interpreter','latex')
set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
title(info.runName)
exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'crests_per_frame_vs_x.pdf'])
%
close all
%
wave_file = [info.rootMat,filesep,info.rootName,'wave_statistics.mat'];
info.waveStatsFile = wave_file;
if exist(wave_file,'file')
    save(wave_file,'-append')
else
     save(wave_file,'-v7.3')
end
save(info.fileName,'-struct','info')
%
