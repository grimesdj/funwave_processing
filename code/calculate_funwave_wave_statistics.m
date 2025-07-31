function [Hs_xy,eta_bar,mask0,x,y,h,fq,Snn_xy,xsl] = calculate_funwave_wave_statistics(rootMat,rootName,bathyFile,Hs,Tp,subDomain,x,y,h);
%
% USAGE: [Hs_x,x,eta_bar,h_bar,sig_eta,freqs,Snn_xy,Snn_x,Snn_wg,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(rootMat,rootName,bathyFile,Hs,Tp,subDomain,plotter,iwg);
%
% plot sea-surface elevation
% need to add some way of figuring out the shoreline location
% xsl = 498;% x=0 shoreline location
%
%
% file_mat = [root,bath,'/',run,'/','dep.mat'];
% h = dep;

% load grid and bottom
if ~exist('h','var')
    load(bathyFile);
end
%
if exist('subDomain','var')
    x = x(subDomain(3):subDomain(4));
    y = y(subDomain(1):subDomain(2));
    h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
end
%
% map bottom points to eta points
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
           h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
nx = length(x);
ny = length(y);
%
% load data
files0 = dir([rootMat,rootName,'mask*.mat']);
files  = dir([rootMat,rootName,'eta*.mat']);
files1 = dir([rootMat,rootName,'u*.mat']);
files2 = dir([rootMat,rootName,'v*.mat']);
if isempty(files2), oneDrun = 1; else oneDrun=0; end
mask0  = [];
eta0   = [];
u0     = [];
v0     = [];
t0     = [];
for ii=1:length(files)
    fprintf('loading eta from: %s \n', files(ii).name);
    if exist('subDomain','var')
        dat0= matfile([files(ii).folder,'/',files0(ii).name]);
        dat = matfile([files(ii).folder,'/',files(ii).name]);
        dat1= matfile([files(ii).folder,'/',files1(ii).name]);
        mask= dat0.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        t   = dat.t;
        eta = dat.eta  (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        u   = dat1.u   (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        if ~oneDrun
            dat2= matfile([files(ii).folder,'/',files2(ii).name]);
            v   = dat2.v   (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        end
    else
        load([files(ii).folder,'/',files0(ii).name]);        
        load([files(ii).folder,'/',files(ii).name]);
        load([files(ii).folder,'/',files1(ii).name]);
        if ~oneDrun
            load([files(ii).folder,'/',files2(ii).name]);
        end
    end
    mask0 = cat(3,mask0,mask);    
    eta0  = cat(3,eta0 ,eta);
    u0    = cat(3,u0   ,u);
    t0    = cat(1,t0   ,t);
    if ~oneDrun
        v0    = cat(3,v0   ,v);
    end
end
eta  = eta0;
u    = u0;
if ~oneDrun
    v= v0;
else
    v= 0*u;
end
t    = t0;
nt   = length(t);
dt   = mean(diff(t));
mask = logical(mask0);
%
% correct for strange numerical behavior in wet-dry region where \eta<h
% eta is referenced to z=<eta>=z_tide; but in swash zone where z_b>=0 the mean water depth depends on eta-z_b,
% which should never be <0;
z_ref  =0-h.*(h<mean(eta(:)));
z_ref_t= repmat(z_ref,1,1,nt);
ht     = repmat(h,1,1,nt);
eta(eta>ht & eta<z_ref_t) = z_ref_t(eta>ht & eta<z_ref_t);
%
clear mask0 eta0 t0 u0 v0
%
% $$$ figure, imagesc(t,x,reshape(eta(1,:,:),nx,nt)),colormap(cmocean('balance')),caxis([-1.5 1.5]),set(gca,'xlim',[500 550],'ylim',x(1)+[0 150]),title('raw $\eta(t,x)$','interpreter','latex')
%
% determine a mean shoreline location
% 5x5 filter
NFx = 5;
NFy = min(5,ny-2); if (NFy<=0), NFy=1; elseif ~mod(NFy,2), NyF=NFy+1; end
filter = ones(NFy,NFx);filter = filter./sum(filter(:));
mask0  = mean(mask,3);
mask0  = conv2([ones(2*NFy+ny,NFx),[ones(NFy,nx);mask0;ones(NFy,nx)],ones(2*NFy+ny,NFx)],filter,'same');
mask0(mask0<1)=0;
mask0 = mask0(NFy+1:end-NFy,NFx+1:end-NFx);
%
% for estimating eta_bar, nan eta
eta(~mask)   = nan;
% $$$ u  (~mask)   = nan;
% $$$ v  (~mask)   = nan;
%
% time-averaged water-level
eta_bar=nanmean(eta,3);
% $$$ figure, imagesc(t,x,reshape(eta(1,:,:),nx,nt)-z_ref(1,:)','AlphaData',reshape(mask(1,:,:),nx,nt)),colormap(cmocean('balance')),caxis([-1.5 1.5]),set(gca,'xlim',[500 550],'ylim',x(1)+[0 150]),title('masked $\eta(t,x)-z_\mathrm{ref}$','interpreter','latex')
%
% now remove mean get Hs across all frequencies
eta        = eta-repmat(eta_bar,1,1,nt);
sig_eta    = 4*nanstd(eta,[],3);
% $$$ figure, imagesc(t,x,reshape(eta(1,:,:),nx,nt),'AlphaData',reshape(mask(1,:,:),nx,nt)),colormap(cmocean('balance')),caxis([-1.5 1.5]),set(gca,'xlim',[500 550],'ylim',x(1)+[0 150]),title('masked $\eta(t,x)-\langle\eta\rangle-z_\mathrm{ref}$','interpreter','latex')
%
% now relax masked region to zero so as to compute spectra in the occationally wetted swash-zone
eta(~mask) = 0;
u          = u.*mask;
v          = v.*mask;
%
% estimate the shoreline location (west coast/east coast)
dum = ones(ny,1)*[1:nx];
if h(1,1)<h(1,end)
    dum((eta_bar)-h-0.12<0)=-inf;
    ixsl = max(dum,[],2);
else
    %    ishoreline = find(h_bar+eta_bar>=0.1,1,'first');
    dum((eta_bar.*mask0)-h-0.12>0)=-inf;
    ixsl = min(dum,[],2);
end
% if there is no shoreline, reference to x=0.
ixsl(isinf(ixsl)) = find(x>=0,1,'first');
xsl = x(ixsl);
%
% compute wave statistics
% with target spectral resolution 1/(n*dt)
df0   = 0.0125;
chnk  = floor(nt*mean(dt)*df0);
inds  = 1:chnk*floor(1/(dt*df0));
%
% reshape to have (t,[x,y])
dum = permute(eta,[3 1 2]);
dum2= reshape(dum, [nt, nx*ny]);
eta_t_xy = dum2(inds,:);
%
dum = permute(u,[3 1 2]);
dum2= reshape(dum, [nt, nx*ny]);
u_t_xy = dum2(inds,:);
%
dum = permute(v,[3 1 2]);
dum2= reshape(dum, [nt, nx*ny]);
v_t_xy = dum2(inds,:);
%
% make depth array to match u,v,eta
h_t_xy = reshape(h,[1,nx*ny]);
%
[Spp,fq]    = mywelch(eta_t_xy,dt,chnk,0.0);
[Suu,~]     = mywelch(u_t_xy  ,dt,chnk,0.0);
[Svv,~]     = mywelch(v_t_xy  ,dt,chnk,0.0);
[Suv,~]     = mycowelch(u_t_xy  ,v_t_xy,dt,chnk,0.0);        
[Spu,~]     = mycowelch(eta_t_xy,u_t_xy,dt,chnk,0.0);
[Spv,~]     = mycowelch(eta_t_xy,v_t_xy,dt,chnk,0.0);
%
df      = fq(2)-fq(1);
fnquist = 0.5/dt;
fq  = fq(2:end);
%
% some nans due to spetra of masked regions make vars complex
SePP= real(Spp(2:end,:));
Suu = real(Suu(2:end,:));
Svv = real(Svv(2:end,:));
Suv = real(Suv(2:end,:));
Spu = real(Spu(2:end,:));
Spv = real(Spv(2:end,:));
%
% convert (u,v) to surface elevation spectra
alpha  = 0.531;
om     = 2*pi*fq;
k      = wavenumber(om,h_t_xy);
cU2eta = (1./om).*real(sinh(k.*h_t_xy)./cosh(k.*h_t_xy.*(1-alpha)));
%
SeUU = Suu.*cU2eta.^2;
SeVV = Svv.*cU2eta.^2;
SeUV = Suv.*cU2eta.^2;
SePU = Spu.*cU2eta;
SePV = Spv.*cU2eta;
%
% estimate bulk stats
% I think there is a sign issue here
coSue =   SePU;
coSve =   SePV;
coUVe =   SeUV;
r2d   = 180/pi;
%
% get a1
a1      = coSue ./ sqrt( SePP .* ( SeUU + SeVV ) );
b1      = coSve ./ sqrt( SePP .* ( SeUU + SeVV ) );
dir1    = r2d* ( atan2(b1,a1) );          
spread1 = r2d* ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );
%
% average over wind-wave band
df       = fq(2)-fq(1);
I        = find(fq>=1/20 & fq<=1/4);
m0       = nansum(SePP(I,:)*df);
m1       = nansum(fq(I).*SePP(I,:)*df);        
ma1      = nansum(a1(I,:).*SePP(I,:)*df,1)/m0;
mb1      = nansum(b1(I,:).*SePP(I,:)*df,1)/m0;
mdir1    = r2d*atan2(mb1,ma1);
mspread1 = r2d*sqrt(abs(2*(1-(ma1.*cos(mdir1/r2d) + mb1.*sin(mdir1/r2d)))));
%
% get a2 b2
a2      = (SePP - SeVV) ./ (SePP + SeVV);
b2      = 2 .* coUVe ./ ( SePP + SeVV );
spread2 = r2d*sqrt(abs(0.5-0.5*(a2.*cos(2.*dir1/r2d)+b2.*sin(2.*dir1/r2d))));
%
ma2      = nansum(a2(I,:).*SePP(I,:)*df,1)/m0;
mb2      = nansum(b2(I,:).*SePP(I,:)*df,1)/m0;
mdir2    = r2d/2*atan2(mb2,ma2);
mspread2 = r2d*sqrt(abs(0.5-0.5*(ma2.*cos(2.*mdir1/r2d)+mb2.*sin(2.*mdir1/r2d))));
%
fq    = fq(I);
Snn_xy= reshape(SePP(I,:),[length(I),ny,nx]);
Hs_xy = reshape(4*sqrt(nansum(df*SePP(I,:))),[ny,nx]);
Tm    = m0./m1;
%
Hs_x   = nanmean(Hs_xy,1);
%
end
