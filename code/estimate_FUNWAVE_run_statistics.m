function info = estimate_FUNWAVE_run_statistics(info);
%
% usage: info = estimate_FUNWAVE_run_statistics(info);
%
% calculate run fast/slow time statistics:
% 1) load (eta,mask,nubrk):
%    1.1) wave frequency spectra, wave height,
%    1.2) breking front statistics
% 2) load (u,v,p,q):
%    2.1) estimate rotational decomposition (Urot,Vrot), and vorticity
%    2.2) estimate breaking dissipation rate and wave force
%    2.3) estimate budget terms:
%    2.3.1) Rotational power: <u_rot dot Fbr>, <u_rot> dot <Fbr>, <u_rot' dot Fbr'>
%                             coh(u_rot,Fbr), 
%    2.3.1) Rotational impulse: <vort dot Fbr>, <vort> dot <Fbr>, <vort' dot Fbr'>
%                              coh(vort,Fbr), coh( int(vort), int(Fbr) )
%
% 3) wave average (30-seconds) fields:
%    (Urot,Vrot, VORT), Fbr, eta 

%    
depFile = [info.rootMat,info.rootName,'dep.nc'];
x0 = ncread(depFile,'x');
y0 = ncread(depFile,'y');
h0 = ncread(depFile,'dep');
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
x = x0(iX);
y = y0(1:ny);
h = h0(subDomain(1):subDomain(2), subDomain(3):subDomain(4));
%
% pre-define file names
info.waveForceFile       = [info.rootMat,info.rootName,'wave_forcing.nc'];
info.waveForceDecompFile = [info.rootMat,info.rootName,'wave_forcing_decomposition.nc'];
%
% load data
files  = dir([info.rootMat,info.rootName,'eta*.nc']);
eta0   = [];
t0     = [];
for ii=1:length(files)
    fprintf('loading eta from: %s \n', files(ii).name);
    fin = [rootMat,rootName,'eta*.nc'];
    eta = ncread(fin,'eta',subDomain([1 3]),subDomain([2 4]));
    t   = ncread(fin,'t');
    %
    dat = ncread(fin,'eta',subDomain([1 3]),subDomain([2 4]));
    % 1) estimate wave forcing, power, vorticity

    % 2) wave average statistics

    % 3) helmholtz decomposition

    % 4) archive
    
    % save eta and time for wave stats calculations
    eta0 = cat(3,eta0,eta);
    t0   = cat(1,t0,t);
end
eta = eta0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
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
