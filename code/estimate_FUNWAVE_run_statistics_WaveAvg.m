function info = estimate_FUNWAVE_run_statistics_WaveAvg(info);
%
% usage: info = estimate_FUNWAVE_run_statistics_WaveAvg(info);
%
% calculate run fast/slow time statistics:
% 1) load (eta,mask,BrkSrc): save these in structure... keep until end!
%    1.1) wave frequency spectra, wave height,
%    1.2) breking front statistics
% 2) load (u_wavg,v_wavg): 
%    2.1) estimate rotational decomposition (Urot,Vrot), and vorticity
%
% 3) Archive wave average (30-seconds) fields: 
%    (Urot,Vrot, VORT), Fbr, eta, and above statistics and archive 

% load the depth and ancillary fields
depFile = [info.rootMat,info.rootName,'dep.nc'];
t0 = ncread(depFile,'t'); nt0= length(t0);
x0 = ncread(depFile,'x'); 
y0 = ncread(depFile,'y'); 
h0 = ncread(depFile,'dep');
dt = t0(2)-t0(1);
dx = x0(2)-x0(1);

% small search radius for front extraction
r0 = floor(2.5/dx);% this is 2.5m in x, and 5m in y for (dx=0.5,dy=1) meters
dy = y0(2)-y0(1);
[xx,yy] = meshgrid(x0,y0);
%
% full alongshore domain
if ~isfield(info,'subDomain')
    iX = find(x0>=0  & x0<=400);
    sdx = length(iX);
    sdy = length(y0);
    subDomain = [1 sdy iX(1) iX(end)];
else
    subDomain = info.subDomain;
end
%
x = x0(subDomain(3):subDomain(4));
y = y0(subDomain(1):subDomain(2));
h = h0(subDomain(1):subDomain(2), subDomain(3):subDomain(4));
nx = length(x);
ny = length(y);
%
% pre-define file names
info.waveForceFile       = [info.rootMat,info.rootName,'wave_forcing.nc'];
info.rotVelFile          = [info.rootMat,info.rootName,'velocity_decomposition.nc'];
info.waveStatsFile       = [info.rootMat,info.rootName,'wave_statistics.nc'];
%
%
% get sea-surface files:
etafiles        = dir([info.rootMat,info.rootName,'eta*.nc']);
eta_wavg_files  = dir([info.rootMat,info.rootName,'etawavg*.nc']);
eta_mean_files  = dir([info.rootMat,info.rootName,'etamean*.nc']);
% remove wavg and mean from instantaneous eta files
etafiles(ismember({etafiles.name}',cat(1,{eta_wavg_files.name}',{eta_mean_files.name}'))) = [];
%
% $$$ % check for time-mean sea-surface... otherwise, use still depth
% $$$ if ~isempty(eta_mean_files)
% $$$     % load the mean sea-surface elevation for depth averaging:    
% $$$     Nf     = length(eta_mean_files);
% $$$     eta_avg= 0;
% $$$     %
% $$$     for ii=1:Nf
% $$$         fprintf('loading eta_mean from: %s \n', etafiles(ii).name);
% $$$         fin = sprintf([info.rootMat,info.rootName,'etamean_%02d.nc'],ii);
% $$$         eta = ncread(fin,'etamean');
% $$$         eta_avg = eta_avg + eta;
% $$$     end
% $$$     eta_avg = eta_avg/Nf;
% $$$     H_avg   = max(h+eta_avg,0.1);
% $$$     clear eta eta_avg
% $$$ else
% $$$     H_avg   = max(h,0.1);
% $$$ end
% $$$ %
% now loop through instantaneous files for stats...
%
Nf     = length(etafiles);
eta0   = [];
t0     = [];
rms_cFbr0 = 0;
%
iter = 1;
for ii=1:Nf
    fprintf('loading eta from: %s \n', etafiles(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'eta_%02d.nc'],ii);
    ncid= netcdf.open(fin);
    [~,nt]  = netcdf.inqDim(ncid,2);
    netcdf.close(ncid)
    eta     = ncread(fin,'eta',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);
    t       = ncread(fin,'t');
    %
    % 
    % 2) estimate wave forcing, power, vorticity
    % 2.0) read in (p,q,brksrc), estimate Fbr,
    vars = {'BrkSrcX','BrkSrcY'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
    end
    H = h+eta;
    % 2.2) estimate breaking dissipation rate and wave force
    %    [Fbx, Fby] = estimate_Fbr(p,q,nubrk,H,dx,dy);
    Fbx      = BrkSrcX./H;
    Fby      = BrkSrcY./H;
    brk_mask = (BrkSrcX ~=0 | BrkSrcY ~=0);
    clear BrkSrc*
    %
    % estimate curl of Fbr
    [Fbx_y,  ~  ] = gradientDG(Fbx./dy);
    [~    ,Fby_x] = gradientDG(Fby./dx);        
    cFbr  = Fby_x - Fbx_y;
    clear Fbx_y Fby_x Fbx Fby
    %
    rms_cFbr = rms(cFbr, 3);
    rms_cFbr0= rms_cFbr0+rms_cFbr.^2;
    %
    eta0  = cat(3,eta0  ,eta);
    clear eta
    % process breaking mask
    for kk = 1:length(t);
        nu    = brk_mask(:,:,kk);
        rclog = bore_front_search_funwave_v2(nu',ny,nx,r0,0.5);
        if isempty(rclog)
            continue
	end
        for ww = 1:length(rclog)
            % convert row/col to x/y
            rc = rclog{ww};
            if isempty(rc)
                continue
            end
            cf = rc(:,2);
	    rf = rc(:,1);
            xf = x(rf);
            yf = y(cf);
            % make sure points are oriented continuously south to north
            [Y,srt]=sort(yf);
            X = xf(srt);
            % estimate front length
% $$$             dl = sqrt(diff(X).^2 + diff(Y).^2);
% $$$             l(iter)  = sum(dl);
            l(iter)  = max(Y)-min(Y);
            xl(iter) = mean(X);
            yl(iter) = mean(Y);
            nl(iter) = length(X);
            tl(iter) = t(kk);
            iter = iter+1;
            xylog{ww,kk} = [X Y];
        end
    end
    %
    % 4) archive eta
    t0   = cat(1,t0,t);
    clear brk_mask t
end
%
rms_cFbr0 = sqrt(rms_cFbr0/Nf);
%
% calculate wave stats... from: estimate_FUNWAVE_run_wave_stats.m
eta = eta0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
clear t0 eta0 nubrk0
%
[Hs_xy,eta_bar,mask0,freqs,Snn_xy,xsl] = estimate_wave_stats_from_eta_h_t_dt(eta,h,t,dt,x);
%
% log shoreline info
info.x_shoreline = xsl;
info.mask        = mask0;
info.eta_bar     = eta_bar;
save(info.fileName,'-struct','info')
%
%
db = 2;
Xbins = [25:db:250];
% keep crest-length stats 
nframes = length(t);
N       = size(xylog,1)/nframes;
mL      = exp(mean(log(l)));
sL      = exp(mean(log(l))+std(log(l)));
Lbins   = [0:10:range(y/2)];
% histogram of lengths
pL      = hist(l,Lbins);
% $$$ mean_stats = struct('N',N,'Length',mL,'Length_plus_std',sL,'Lbins',Lbins,'Length_histogram',pL);
%
% cross-shore bins
Nx  = nan*Xbins;
mLx = nan*Xbins;
sLx = nan*Xbins;
Lx_log_mean = nan*Xbins;
Lx_log_std  = nan*Xbins;
Lx_histo    = nan*(Lbins'*Xbins);
for bin = 1:length(Xbins)
    iX = find(xl>=Xbins(bin)-db/2 & xl<Xbins(bin)+db/2);
    Nx(bin) = length(iX)/nframes;
    Lx_log_mean(bin)=mean(log(l(iX)));
    Lx_log_std(bin)=std(log(l(iX)));
    mLx(bin)=exp(mean(log(l(iX))));
    sLx(bin)=exp(log(mLx(bin))+std(log(l(iX))));
    tmp = hist(l(iX),Lbins);
    Lx_histo(:,bin) = tmp;
end
%
%
% binned_stats = struct('Xbins',Xbins,'N',Nx,'Length',mLx,'Length_plus_std',sLx,'Length_histogram',Lx_histo,'log_mean_length',Lx_log_mean,'log_std_length',Lx_log_std);
%
%
%% Calculate rotational current decomposition
files  = dir([info.rootMat,info.rootName,'uwavg*.nc']);
Nf     = length(files);
t0     = [];
Urot   = [];
Vrot   = [];
PSI    = [];
VORT   = [];
ERR    = [];
eta0   = [];
iter = 1;
for ii=1:Nf
    %
    fprintf('loading eta from: %s \n', files(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'uwavg_%02d.nc'],ii);
    t   = ncread(fin,'t');
    nt  = length(t);
    %
    % 3) vorticity statistics
    vars = {'etawavg','uwavg','vwavg'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
    end
    %
    eta=etawavg;
    u = uwavg;
    v = vwavg;
    clear etawavg uwavg vwavg
    %
    % estimate vorticity
    [uy,~ ,~] = gradientDG(u./dy);
    [~ ,vx,~] = gradientDG(v./dx);
    omega = vx-uy;
    clear uy vx
    %
    %
    %
    % helmholtz decomposition
    II=0;
    for jj = 1:nt
        II=II+1;
        U = u(:,:,jj);
        V = v(:,:,jj);
        [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(U,V,dx,dy);
        urot(:,:,II) = u_psi;
        vrot(:,:,II) = v_psi;
        psi0(:,:,II) = psi;
        err (:,:,II) = sqrt( (U-(u_psi+u_phi)).^2 + (V-(v_psi+v_phi)).^2 );
    end
    % log Urot, Vrot
    t0   = cat(1,t0,t);
    VORT = cat(3,VORT,omega);
    Urot = cat(3,Urot,urot);
    Vrot = cat(3,Vrot,vrot);
    PSI = cat(3,PSI,psi0);
    ERR = cat(3,ERR,err);
    eta0= cat(3,eta0,eta);
    clear err psi0 urot vrot u_psi v_psi U V 
    %
    %
end
%
%
save(info.fileName,'-struct','info')
%
% For debugging the netcdf write portion, uncomment this:
% $$$ save('/scratch/grimesdj/ripchannel/planar2D/mat_data/debugging_code.mat','-v7.3')
%
% Archive mean fields to netcdf files:
% 1) info.waveForceFile       = [info.rootMat,info.rootName,'wave_forcing.nc'];
% 1.1) variables:
%      Eb0   = high-pass dissipation(?)  ( u' \dot Fbr )
%      Ebavg =  low-pass generation(?)   (\avg{u} \dot Fbr)
%      Icoh   =  coherence between diff(vort) and int(Fbr*dt)
%      Icoh_avg =  coherence between wave averaged diff(vort) and int(Fbr*dt)
%      Xbins =  cross-shore bins for crest-length statistics
%      Nc    =  number of crests per bin (Nx)
%      Lc    =  mean crest lengt per bin (mLx)
%      Lc_std=  mean + std per bin (sLx)
%      Lbins =  bins for crest length PDF
%      Lc_pdf=  histogram of lengths (Lx_histo)
if exist(info.waveForceFile,'file')
    eval(['!rm ',info.waveForceFile])
end
% this is for archiving
dim_yx = {"y",length(y),"x",length(x)};
nccreate  (info.waveForceFile,'rms_cFbr','Dimensions',dim_yx,'Format','netcdf4')
ncwrite   (info.waveForceFile,'rms_cFbr',rms_cFbr0);
ncwriteatt(info.waveForceFile,'rms_cFbr','Description','rms of the curl of viscous breaking force');
%
nccreate  (info.waveForceFile,'x','Dimensions',{"x",length(x)},'Format','netcdf4')
ncwrite   (info.waveForceFile,'x',x);
nccreate  (info.waveForceFile,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
ncwrite   (info.waveForceFile,'y',y);
%
dim_lbxb = {"Lb",length(Lbins),"Xb",length(Xbins)};
nccreate  (info.waveForceFile,'Lc_pdf','Dimensions',dim_lbxb,'Format','netcdf4')
ncwrite   (info.waveForceFile,'Lc_pdf',Lx_histo./sum(Lx_histo,1));
ncwriteatt(info.waveForceFile,'Lc_pdf','Description','cross-shore bin averaged crest-length pdf');
%
nccreate  (info.waveForceFile,'Xb','Dimensions',{"Xb",length(Xbins)},'Format','netcdf4')
ncwrite   (info.waveForceFile,'Xb',Xbins);
nccreate  (info.waveForceFile,'Lb','Dimensions',{"Lb",length(Lbins)},'Format','netcdf4')
ncwrite   (info.waveForceFile,'Lb',Lbins);
%
dim_xb  = {"Xb",length(Xbins)};
nccreate  (info.waveForceFile,'Lc','Dimensions',dim_xb,'Format','netcdf4')
ncwrite   (info.waveForceFile,'Lc',mLx);
ncwriteatt(info.waveForceFile,'Lc','Description','cross-shore bin averaged crest length');
%
nccreate  (info.waveForceFile,'Nc','Dimensions',dim_xb,'Format','netcdf4')
ncwrite   (info.waveForceFile,'Nc',Nx);
ncwriteatt(info.waveForceFile,'Nc','Description','cross-shore bin averaged number of crests');
%
nccreate  (info.waveForceFile,'Lc_plus_std','Dimensions',dim_xb,'Format','netcdf4')
ncwrite   (info.waveForceFile,'Lc_plus_std',sLx);
ncwriteatt(info.waveForceFile,'Lc_plus_std','Description','cross-shore bin averaged crest length plus one standard deviation');
%
% 2) info.rotVelFile          = [info.rootMat,info.rootName,'velocity_decomposition.nc'];
% 2.1) variables:
%      VORT = low-pass vorticity
%      Urot = low-pass rotational cross-shore velocity
%      Vrot = low-pass rotational alongshore velocity
%      PSI  = low-pass streamfunction
%      ERR  = error between decomposition and actual
if exist(info.rotVelFile,'file')
    eval(['!rm ',info.rotVelFile])
end
dim_yxt = {"y",length(y),"x",length(x),"t",length(t0)};
nccreate  (info.rotVelFile,'VORT','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'VORT',VORT);
ncwriteatt(info.rotVelFile,'VORT','Description','low-pass vorticity');
%
nccreate  (info.rotVelFile,'x','Dimensions',{"x",length(x)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'x',x);
nccreate  (info.rotVelFile,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'y',y);
nccreate  (info.rotVelFile,'t','Dimensions',{"t",length(t0)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'t',t0);
%
nccreate  (info.rotVelFile,'Urot','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'Urot',Urot);
ncwriteatt(info.rotVelFile,'Urot','Description','low-pass rotational cross-shore velocity');
%
nccreate  (info.rotVelFile,'Vrot','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'Vrot',Vrot);
ncwriteatt(info.rotVelFile,'Vrot','Description','low-pass rotational alongshore velocity');
%
nccreate  (info.rotVelFile,'PSI','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'PSI',PSI);
ncwriteatt(info.rotVelFile,'PSI','Description','low-pass streamfunction');
%
nccreate  (info.rotVelFile,'ERR','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'ERR',ERR);
ncwriteatt(info.rotVelFile,'ERR','Description','low-pass error in helmholtz decomposition');
%
nccreate  (info.rotVelFile,'eta','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'eta',eta0);
ncwriteatt(info.rotVelFile,'eta','Description','low-pass sea-surface elevation');
%
% 3) info.waveStatsFile       = [info.rootMat,info.rootName,'wave_statistics.nc'];
% 3.1) variables:
%     Hs_xy  = spacially varying significant wave height
%     eta_bar= time averaged waterlevel
%     mask   = mask based on time-averaged waterlevel
%     xsl    = shoreline location based on mask
%     freqs  = frequencies 
%     Snn_xy = waterlevel spectrum versus space
if exist(info.waveStatsFile,'file')
    eval(['!rm ',info.waveStatsFile])
end
dim_yx  = {"y",length(y),"x",length(x)};
nccreate  (info.waveStatsFile,'Hs','Dimensions',dim_yx,'Format','netcdf4')
ncwrite   (info.waveStatsFile,'Hs',Hs_xy);
ncwriteatt(info.waveStatsFile,'Hs','Description','Wave height');
%
nccreate  (info.waveStatsFile,'x','Dimensions',{"x",length(x)},'Format','netcdf4')
ncwrite   (info.waveStatsFile,'x',x);
nccreate  (info.waveStatsFile,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
ncwrite   (info.waveStatsFile,'y',y);
%
nccreate  (info.waveStatsFile,'eta','Dimensions',dim_yx,'Format','netcdf4')
ncwrite   (info.waveStatsFile,'eta',eta_bar);
ncwriteatt(info.waveStatsFile,'eta','Description','Time-averaged waterlevel');
%
nccreate  (info.waveStatsFile,'mask','Dimensions',dim_yx,'Format','netcdf4')
ncwrite   (info.waveStatsFile,'mask',mask0);
ncwriteatt(info.waveStatsFile,'mask','Description','land mask based on time averaged waterlevel');
%
dim_y  = {"y",length(y)};
nccreate  (info.waveStatsFile,'xsl','Dimensions',dim_y,'Format','netcdf4')
ncwrite   (info.waveStatsFile,'xsl',xsl);
ncwriteatt(info.waveStatsFile,'xsl','Description','shoreline coordinates based on time averaged waterlevel');
%
dim_fx = {"frq",length(freqs),"y",length(y),"x",length(x)};
nccreate  (info.waveStatsFile,'Snn','Dimensions',dim_fx,'Format','netcdf4')
ncwrite   (info.waveStatsFile,'Snn',Snn_xy);
ncwriteatt(info.waveStatsFile,'Snn','Description','waterlevel spectra versus along- and cross-shore coordinate');
%
nccreate  (info.waveStatsFile,'frq','Dimensions',{"frq",length(freqs)},'Format','netcdf4')
ncwrite   (info.waveStatsFile,'frq',freqs);
%
return
% $$$ 
% $$$ 
% $$$ 
% $$$ % this if from the wave stats code
% $$$ eta = eta0;
% $$$ t   = t0;
% $$$ nt  = length(t);
% $$$ dt  = mean(diff(t));
% $$$ %
% $$$ [Hs_xy,eta_bar,mask0,x,y,h,freq,Snn_xy,xsl] = calculate_funwave_wave_height_statistics_v2(info.rootMat,info.rootName,info.bathyFile,info.Hs,info.Tp,subDomain);
% $$$ %
% $$$ info.x_shoreline = xsl;
% $$$ info.mask        = mask0;
% $$$ info.eta_bar     = eta_bar;
% $$$ %
% $$$ fig0 = figure;
% $$$ Hs_x = nanmean(Hs_xy,1);
% $$$ p0 = plot(x,Hs_xy,'.k',x,Hs_x,'-r','markersize',1,'linewidth',1.5);
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$H_s$ [m]','interpreter','latex')
% $$$ set(gca,'xlim',[75 600],'ylim',[0 1.1*nanmax(Hs_xy(:))],'ticklabelinterpreter','latex','tickdir','out')
% $$$ f0l1 = legend([p0(1),p0(end)]','$H_\mathrm{s}(x,y)$','$\bar{H}_\mathrm{s}(x)$');
% $$$ set(f0l1,'location','southeast','interpreter','latex')
% $$$ title(info.runName)
% $$$ if ~exist([info.rootSim,filesep,'figures'],'dir')
% $$$     eval(['!mkdir ',[info.rootSim,filesep,'figures']])
% $$$ end
% $$$ exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'Hs.pdf'])
% $$$ %
% $$$ end
% $$$ %
% $$$ % get front statistics
% $$$ [mean_stats,binned_stats] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ % $$$ [mL,sL,N,mLx,sLx,Lx_log_mean,Lx_log_std,Nx,xylog,Xbins] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ % front_file = [info.rootMat,filesep,info.rootName,'bore_front_statistics.mat'];
% $$$ %
% $$$ % make a stats plot? what stats?
% $$$ fig1 = figure;
% $$$ p1 = plot(binned_stats.Xbins,binned_stats.Length,'-k',binned_stats.Xbins,exp( binned_stats.log_mean_length' + binned_stats.log_std_length'*[-1 1] ), '--r');
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$L$ [m]','interpreter','latex')
% $$$ set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
% $$$ f1l1 = legend([p1(1), p1(2)]','$\bar{L}(l)=\exp\left(\overline{\log(l)}\right)$','$\bar{L}\pm\mathrm{std}(L)$');
% $$$ set(f1l1,'location','northeast','interpreter','latex')
% $$$ title(info.runName)
% $$$ if ~exist([info.rootSim,filesep,'figures'],'dir')
% $$$     eval(['!mkdir ',[info.rootSim,filesep,'figures']])
% $$$ end
% $$$ exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'crest_length_vs_x.pdf'])
% $$$ %
% $$$ fig2 = figure;
% $$$ p2 = plot(binned_stats.Xbins,binned_stats.N,'-k','linewidth',2);
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$N$ [crests/frame]','interpreter','latex')
% $$$ set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
% $$$ title(info.runName)
% $$$ exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'crests_per_frame_vs_x.pdf'])
% $$$ %
% $$$ close all
% $$$ %
% $$$ wave_file = [info.rootMat,filesep,info.rootName,'wave_statistics.mat'];
% $$$ info.waveStatsFile = wave_file;
% $$$ if exist(wave_file,'file')
% $$$     save(wave_file,'-append')
% $$$ else
% $$$      save(wave_file,'-v7.3')
% $$$ end
% $$$ save(info.fileName,'-struct','info')
% $$$ %
