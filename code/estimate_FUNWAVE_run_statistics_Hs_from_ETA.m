function info = estimate_FUNWAVE_run_statistics_Hs_from_ETA(info);
%
% usage: info = estimate_FUNWAVE_run_statitics_Waves(info);
%
% calculate run fast-time wave statistics:
% 1) load (eta,mask,BrkSrc): save these in structure... keep until end!
%    1.1) wave frequency spectra, wave height,
%    1.2) breking front statistics

% load the depth and ancillary fields
depFile = [info.rootMat,info.rootName,'dep.nc'];
t0 = ncread(depFile,'t'); nt0= length(t0);
x0 = ncread(depFile,'x'); 
y0 = ncread(depFile,'y'); 
h0 = ncread(depFile,'dep');
dt = t0(2)-t0(1);
dx = x0(2)-x0(1);

% $$$ % small search radius for front extraction
% $$$ r0 = floor(2.5/dx);% this is 2.5m in x, and 5m in y for (dx=0.5,dy=1) meters
% $$$ dy = y0(2)-y0(1);
% $$$ [xx,yy] = meshgrid(x0,y0);
% $$$ %
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
% $$$ info.waveForceFile       = [info.rootMat,info.rootName,'wave_forcing.nc'];
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
% load mean sealevel from momentum file:
%
% now loop through instantaneous files for stats...
%
Nf     = length(etafiles);
eta1   = 0;
eta2   = 0;
% $$$ eta0   = [];
% $$$ t0     = [];
% $$$ rms_cFbr0 = 0;
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
    nt      = length(t);
    H       = h0+eta;
    eta(H<0.1)=0;
    eta2 = eta2+mean(eta.^2,3);
    eta1 = eta1+mean(eta   ,3);
    %
    % 
% $$$     % 2) estimate wave forcing, power, vorticity
% $$$     % 2.0) read in (p,q,brksrc), estimate Fbr,
% $$$     vars = {'BrkSrcX','BrkSrcY'};
% $$$     for jj=1:length(vars)
% $$$         fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
% $$$         eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
% $$$     end
% $$$     H = h+eta;
% $$$     % 2.2) estimate breaking dissipation rate and wave force
% $$$     %    [Fbx, Fby] = estimate_Fbr(p,q,nubrk,H,dx,dy);
% $$$     Fbx      = BrkSrcX./H;
% $$$     Fby      = BrkSrcY./H;
% $$$     brk_mask = (BrkSrcX ~=0 | BrkSrcY ~=0);
% $$$     clear BrkSrc*
% $$$     %
% $$$     % estimate curl of Fbr
% $$$     [Fbx_y,  ~  ] = gradientDG(Fbx./dy);
% $$$     [~    ,Fby_x] = gradientDG(Fby./dx);        
% $$$     cFbr  = Fby_x - Fbx_y;
% $$$     clear Fbx_y Fby_x Fbx Fby
% $$$     %
% $$$     rms_cFbr = rms(cFbr, 3);
% $$$     rms_cFbr0= rms_cFbr0+rms_cFbr.^2;
% $$$     %
% $$$     eta0  = cat(3,eta0  ,eta);
% $$$     clear eta
% $$$     % process breaking mask
% $$$     for kk = 1:length(t);
% $$$         nu    = brk_mask(:,:,kk);
% $$$         rclog = bore_front_search_funwave_v2(nu',ny,nx,r0,0.5);
% $$$         if isempty(rclog)
% $$$             continue
% $$$ 	end
% $$$         for ww = 1:length(rclog)
% $$$             % convert row/col to x/y
% $$$             rc = rclog{ww};
% $$$             if isempty(rc)
% $$$                 continue
% $$$             end
% $$$             cf = rc(:,2);
% $$$ 	    rf = rc(:,1);
% $$$             xf = x(rf);
% $$$             yf = y(cf);
% $$$             % make sure points are oriented continuously south to north
% $$$             [Y,srt]=sort(yf);
% $$$             X = xf(srt);
% $$$             % estimate front length
% $$$ % $$$             dl = sqrt(diff(X).^2 + diff(Y).^2);
% $$$ % $$$             l(iter)  = sum(dl);
% $$$             l(iter)  = max(Y)-min(Y);
% $$$             xl(iter) = mean(X);
% $$$             yl(iter) = mean(Y);
% $$$             nl(iter) = length(X);
% $$$             tl(iter) = t(kk);
% $$$             iter = iter+1;
% $$$             xylog{ww,kk} = [X Y];
% $$$         end
% $$$     end
% $$$     %
% $$$     % 4) archive eta
% $$$     t0   = cat(1,t0,t);
% $$$     clear brk_mask t
end
eta2 = eta2/Nf;
eta1 = eta1/Nf;
Hsig = 4*sqrt( eta2 - eta1.^2 );
%
% $$$ rms_cFbr0 = sqrt(rms_cFbr0/Nf);
% $$$ %
% $$$ % calculate wave stats... from: estimate_FUNWAVE_run_wave_stats.m
% $$$ eta = eta0;
% $$$ t   = t0;
% $$$ nt  = length(t);
% $$$ dt  = mean(diff(t));
% $$$ clear t0 eta0 nubrk0
% $$$ %
% $$$ [Hs_xy,eta_bar,mask0,freqs,Snn_xy,xsl] = estimate_wave_stats_from_eta_h_t_dt(eta,h,t,dt,x);
% $$$ %
% $$$ % log shoreline info
% $$$ info.x_shoreline = xsl;
% $$$ info.mask        = mask0;
% $$$ info.eta_bar     = eta_bar;
% $$$ save(info.fileName,'-struct','info')
% $$$ %
% $$$ %
% $$$ db = 2;
% $$$ Xbins = [25:db:250];
% $$$ % keep crest-length stats 
% $$$ nframes = length(t);
% $$$ N       = size(xylog,1)/nframes;
% $$$ mL      = exp(mean(log(l)));
% $$$ sL      = exp(mean(log(l))+std(log(l)));
% $$$ Lbins   = [0:10:range(y/2)];
% $$$ % histogram of lengths
% $$$ pL      = hist(l,Lbins);
% $$$ % $$$ mean_stats = struct('N',N,'Length',mL,'Length_plus_std',sL,'Lbins',Lbins,'Length_histogram',pL);
% $$$ %
% $$$ % cross-shore bins
% $$$ Nx  = nan*Xbins;
% $$$ mLx = nan*Xbins;
% $$$ sLx = nan*Xbins;
% $$$ Lx_log_mean = nan*Xbins;
% $$$ Lx_log_std  = nan*Xbins;
% $$$ Lx_histo    = nan*(Lbins'*Xbins);
% $$$ for bin = 1:length(Xbins)
% $$$     iX = find(xl>=Xbins(bin)-db/2 & xl<Xbins(bin)+db/2);
% $$$     Nx(bin) = length(iX)/nframes;
% $$$     Lx_log_mean(bin)=mean(log(l(iX)));
% $$$     Lx_log_std(bin)=std(log(l(iX)));
% $$$     mLx(bin)=exp(mean(log(l(iX))));
% $$$     sLx(bin)=exp(log(mLx(bin))+std(log(l(iX))));
% $$$     tmp = hist(l(iX),Lbins);
% $$$     Lx_histo(:,bin) = tmp;
% $$$ end
% $$$ %
% $$$ %
% $$$ %
% $$$ save(info.fileName,'-struct','info')
% $$$ %
% For debugging the netcdf write portion, uncomment this:
% $$$ save('/scratch/grimesdj/ripchannel/planar2D/mat_data/debugging_code.mat','-v7.3')
%
momFile  = dir([info.rootMat,info.rootName,'MomentumTerms.nc']);
momFile  = [momFile(1).folder, filesep, momFile.name];
dim_yx = {"y",length(y),"x",length(x)};
finfo = ncinfo(momFile);
vars  = {finfo.Variables.Name}';
if ~ismember('Hsig',vars)
    nccreate  (momFile,'Hsig','Dimensions',dim_yx,'Format','netcdf4')
    ncwriteatt(momFile,'Hsig','Description','Significant wave height estimated from: 4*sqrt( \bar{eta^2}-\bar{eta}^2');
    ncwrite   (momFile,'Hsig',Hsig);
else
    try
        ncwrite   (momFile,'Hsig',Hsig);
    catch
        disp(['cannot add Hsig to file: ',momFile])
    end
end


return
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
