function info = estimate_FUNWAVE_run_wave_forcing(info);
%
% usage:  info = estimate_FUNWAVE_run_wave_forcing(info);
%
% calculate wave breaking force and perform rotational/irrotational decomposition.
% Analysis is limited to the region:
%       75 m <= x <=  525 m
%     -200 m <= y <= 1300 m

% get the input (x,y,h) grid file
load(info.bathyFile)
% switch from h=z_bottom to h=depth
h =-h;
% save the original grid
x0=x;
y0=y;
h0=h;
%
% depth threshold
threshold = 0.01;
%
% full alongshore domain
iX = find(x0>=75  & x0<=525);
subDomain = [1 length(y) iX(1) iX(end)];
%
x = x(subDomain(3):subDomain(4));
y = y(subDomain(1):subDomain(2));
h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
%
% map bottom points to eta points
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
           h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
%
% needed for estimating gradients
[xx,yy] = meshgrid(x,y);
dx = x(2)-x(1);
dy = y(2)-y(1);
nx = length(x);
ny = length(y);
%
% loop over files and only archive what you need?
files0 = dir([info.rootMat,info.rootName,'mask*.mat']);
files1  = dir([info.rootMat,info.rootName,'nubrk*.mat']);
files2  = dir([info.rootMat,info.rootName,'p*.mat']);
files3  = dir([info.rootMat,info.rootName,'q*.mat']);
files4  = dir([info.rootMat,info.rootName,'eta*.mat']);
Nf      = length(files0);
%
FORCE   = matfile([info.rootMat,info.rootName,'wave_forcing.mat'],'writable',true);
DECOMP  = matfile([info.rootMat,info.rootName,'wave_forcing_decomposition.mat'],'writable',true);
info.waveForceFile = [info.rootMat,info.rootName,'wave_forcing.mat'];
info.waveForceDecompFile = [info.rootMat,info.rootName,'wave_forcing_decomposition.mat'];
save(info.fileName,'-struct','info')
%
inds = 0;
FORCE.x = x;
FORCE.y = y;
DECOMP.x = x;
DECOMP.y = y;
cFbr = [];
Havg = 0;
Hmin = inf;
for ii=1:Nf
    fprintf('loading nubrk from: %s \n', files1(ii).name);
    if exist('subDomain','var')
% $$$         dat0 = matfile([files0(ii).folder,'/',files0(ii).name]);
        dat1 = matfile([files1(ii).folder,'/',files1(ii).name]);
        dat2 = matfile([files2(ii).folder,'/',files2(ii).name]);
        dat3 = matfile([files3(ii).folder,'/',files3(ii).name]);
        dat4 = matfile([files4(ii).folder,'/',files4(ii).name]);                        
% $$$         mask = dat0.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        nubrk= dat1.nubrk(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        P    = dat2.p    (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        Q    = dat3.q    (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        eta  = dat4.eta  (subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);        
        t    = dat1.t;
    else
% $$$         load([files0(ii).folder,'/',files0(ii).name]);        
        load([files1(ii).folder,'/',files1(ii).name]);
        load([files2(ii).folder,'/',files2(ii).name]);
        load([files3(ii).folder,'/',files3(ii).name]);
        load([files4(ii).folder,'/',files4(ii).name]);        
    end
    %
    nt = length(t);
    dt = nanmean(diff(t));
    % average to eta/u/v-points?
    H             = eta + repmat(h,1,1,nt);
    % gradient of fluxes
    [Py, Px, ~]   = gradientDG(P);% grad(H.*u);
    [Qy, Qx, ~]   = gradientDG(Q);% grad(H.*v);
    % divergence of fluxes... stress
    [Pyy,   ~, ~] = gradientDG(nubrk.*Py/dy);
    [~  , Pxx, ~] = gradientDG(nubrk.*Px/dx);
    [Qyy,   ~, ~] = gradientDG(nubrk.*Qy/dy);
    [~  , Qxx, ~] = gradientDG(nubrk.*Qx/dx);
    %
    clear Py Px Qy Qx
    Pyy           = Pyy/dy; 
    Pxx           = Pxx/dx;
    Qyy           = Qyy/dy;
    Qxx           = Qxx/dx;
    %
    % body forcing / per unit volume
    Rbx             = (Pyy+Pxx)./H; clear Pxx Pyy
    Rby             = (Qxx+Qyy)./H; clear Qxx Qyy
    [Rbx_y,~    ,~] = gradientDG(Rbx);
    [~    ,Rby_x,~] = gradientDG(Rby);
    %
    % curl of body forcing
    curlFbr = (Rby_x/dx - Rbx_y/dy);
    %
    % apply mask
    mask = (H>threshold);
    Rbx(~mask)     = nan;
    Rby(~mask)     = nan;
    curlFbr(~mask) = nan;
    cFbr = cat(3,cFbr,curlFbr);
    Havg = Havg+sum(H,3)/nt;
    Hmin = min(Hmin,min(H,[],3));
    %
    % which indices are we on?
    inds = inds(end) + [1:nt];
    FORCE.t (inds,1) = t;
    FORCE.H (1:ny,1:nx,inds) = H;
    FORCE.Fx(1:ny,1:nx,inds) = Rbx;
    FORCE.Fy(1:ny,1:nx,inds) = Rby;
    FORCE.curlF(1:ny,1:nx,inds) = curlFbr;    
    %
    for jj=1:size(Rbx,3);
        Fx = Rbx(:,:,jj);
        Fy = Rby(:,:,jj);
        %
        [psi,psi_y,psi_x,phi,phi_x,phi_y]=get_vel_decomposition_reGRID(Fx,Fy,dx,dy);
        PSI(:,:,jj) = psi;
        Rx_rot(:,:,jj)= psi_y;
        Ry_rot(:,:,jj)= psi_x;
        PHI(:,:,jj) = phi;
        Rx_irr(:,:,jj)= phi_x;
        Ry_irr(:,:,jj)= phi_y;
        err(:,:,jj)   = sqrt( (Fx-(psi_y+phi_x)).^2 + (Fy-(psi_x+phi_y)).^2);
    end
    PSI   (~mask)=nan;
    Rx_rot(~mask)=nan;
    Ry_rot(~mask)=nan;
    PHI   (~mask)=nan;
    Rx_irr(~mask)=nan;
    Ry_irr(~mask)=nan;
    err   (~mask)=nan;
    DECOMP.t  (inds,1)         = t;
    DECOMP.H  (1:ny,1:nx,inds) = H;
    DECOMP.PSI(1:ny,1:nx,inds) = PSI;
    DECOMP.Fx_rot(1:ny,1:nx,inds) = Rx_rot;
    DECOMP.Fy_rot(1:ny,1:nx,inds) = Ry_rot;    
    DECOMP.PHI(1:ny,1:nx,inds) = PHI;    
    DECOMP.Fx_irr(1:ny,1:nx,inds) = Rx_irr;    
    DECOMP.Fy_irr(1:ny,1:nx,inds) = Ry_irr;
    DECOMP.err(1:ny,1:nx,inds)    = err;
    clear PSI PHI Rx_rot Ry_rot Rx_irr Ry_irr err
end
%
Havg     = Havg/Nf;
mask     = info.mask;
% time average and rms
cFbr_tavg= nanmean(cFbr,3).*mask;
rms_cFbr_t = sqrt(nanmean((cFbr-cFbr_tavg).^2,3)).*mask;
% alongshore average using x_shoreline from wave_stats analysis
xsl= info.x_shoreline;
avg_xsl = round(mean(xsl));
avg_cFbr = 0*x;
rms_cFbr = 0*x;
for jj = 1:ny
    tmp1 = interp1((x-xsl(jj)), cFbr_tavg(jj,:), x-avg_xsl);
    tmp2 = interp1((x-xsl(jj)), rms_cFbr_t(jj,:).^2, x-avg_xsl);
    avg_cFbr = avg_cFbr+tmp1;
    rms_cFbr = rms_cFbr+tmp2;
end
avg_cFbr   = avg_cFbr/ny;
rms_cFbr   = sqrt(rms_cFbr/ny);
%
FORCE.mask  = mask;
FORCE.Hmin  = threshold;
DECOMP.mask = mask;
DECOMP.Hmin = threshold;
%
fig0 = figure;
p0 = semilogy(x,abs(avg_cFbr),'-k',x,rms_cFbr,'-r');
xlabel('cross-shore [m]','interpreter','latex')
ylabel('$\nabla\times\vec{F}_\mathrm{br}$','interpreter','latex')
f0l0 = legend(p0,'$|\overline{(\cdot)}|$','rms($\cdot$)');
set(f0l0,'interpreter','latex')
set(gca,'xlim',[75 300],'ticklabelinterpreter','latex')
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_forcing_vs_x.pdf'])
%
cFbr_range = range(cFbr_tavg(:));
clims    = [-0.1 0.1];
clrs     = clims(1):diff(clims)/255:clims(2);
cm       = cmocean('balance');
fig1 = figure;
imagesc(x,y,cFbr_tavg,clims),colormap(cm),caxis(clims)
xlabel('cross-shore [m]','interpreter','latex')
ylabel('alongshore [m]','interpreter','latex')
set(gca,'xlim',[75 500],'ticklabelinterpreter','latex')
%
a3 = axes('position',[0.85 0.7 0.05 0.25]);
imagesc(0,clrs,reshape(cm,256,1,3))
ylabel(a3,'$\nabla\times\vec{F}_\mathrm{br}$','interpreter','latex','horizontalalignment','right')
set(a3,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xtick',[])
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_forcing_time_averaged.pdf'])
close all
% close fig0 fig1
%
% $$$ wave_file = [info.rootMat,filesep,info.rootName,'wave_forcing.mat'];
% $$$ info.waveForceFile = wave_file;
% $$$ save(info.fileName,'-struct','info')

end

