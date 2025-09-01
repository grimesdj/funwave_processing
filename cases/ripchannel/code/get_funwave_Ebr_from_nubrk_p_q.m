function [Ebr] = get_funwave_Ebr_from_nubrk_p_q(runNAME,rawDIR,dx,dy,Navg);
%
% USAGE: [Ebr] = get_funwave_Ebr_from_nubrk_p_q(runNAME,rawDIR,dx,dy,Navg);
%
% runNAME: used for archiving stats
% rawDIR:  location of raw ascii data files
% (dy,dx): grid spacing (meters)
% Ebr:     (u-u_avg)\cdot Fbr
% u_avg:   backward average over Navg points
% Fbr:     1/(h+eta)* div( nubrk*(h+eta)*grad( u ) )


% get structure with all elevation data
fprintf('working on Ebr for run: %s \n',runNAME);
eta_files   = dir([rawDIR,filesep,'eta_*']);
nubrk_files = dir([rawDIR,filesep,'nubrk_*']);
p_files     = dir([rawDIR,filesep,'p_*']);
q_files     = dir([rawDIR,filesep,'q_*']);
Nf          = length(eta_files);
%
% get bathymetry file for filtering by water-depth
h      = load([rawDIR,filesep,'dep.out'],'-ascii');
%h      = 0.25*( h(1:end-1,1:end-1) + h(2:end,1:end-1) +...
%		h(2:end  ,1:end-1) + h(2:end,2:end) );
h1yavg = mean(h,1);
[Ny,Nx]= size(h);
%
% pre-allocate variables
N     = 0;
E1sum = zeros(1,Nx);
E2sum = zeros(1,Nx);
mask  = zeros(1,Nx);
%
% loop over Nf files, load and compute stats
for ii = 1:Nf
    % load the eta_ file
    %    fprintf('\tloading files associated with: %s \n',eta_files(ii).name);
    eta      = load([rawDIR,filesep,eta_files(ii).name],'-ascii');
    nubrk    = load([rawDIR,filesep,nubrk_files(ii).name],'-ascii');
    p        = load([rawDIR,filesep,p_files(ii).name],'-ascii');
    q        = load([rawDIR,filesep,q_files(ii).name],'-ascii');
    %
    % make wave-avg (p,q) for estimate of dissipation rate
    if ii<=Navg
        p_t(:,:,ii) = p;
        q_t(:,:,ii) = q;
        continue
    end
    % water depth
    H        = h+eta;
    H(H<=0)  = 0.001;
    % gradient of fluxes
    [Py, Px]   = gradientDG(p);% grad(H.*u);
    [Qy, Qx]   = gradientDG(q);% grad(H.*v);
    % divergence of fluxes... stress
    [Pyy,   ~] = gradientDG(nubrk.*Py/dy);
    [~  , Pxx] = gradientDG(nubrk.*Px/dx);
    [Qyy,   ~] = gradientDG(nubrk.*Qy/dy);
    [~  , Qxx] = gradientDG(nubrk.*Qx/dx);
    %
    clear Py Px Qy Qx
    Pyy           = Pyy/dy; 
    Pxx           = Pxx/dx;
    Qyy           = Qyy/dy;
    Qxx           = Qxx/dx;
    % body forcing / per unit volume
    Fbx             = (Pyy+Pxx)./H; clear Pxx Pyy
    Fby             = (Qxx+Qyy)./H; clear Qxx Qyy
    % wave dissipation rate (u' \cdot Fb)
    Ebr      = Fbx.*(p-mean(p_t,3,'omitnan')) + Fby.*(q-mean(q_t,3,'omitnan'));
    Ebr1yavg = nanmean(Ebr,1);
    Ebr2yavg = nanmean(Ebr.^2,1);
    E1sum    = E1sum + Ebr1yavg;
    E2sum    = E2sum + Ebr2yavg;
    %    mask     = mask  + double(eta1yavg+h1yavg>=0.1);
    N        = N     + 1;
    %
    % cycle the time history and over-write oldest field
    p_t      = circshift(p_t,[0 0 1]); p_t(:,:,Navg) = p;
    q_t      = circshift(q_t,[0 0 1]); q_t(:,:,Navg) = q;
end
sig2= Ny/( (N-1)*(Ny-1)) * (E2sum - 2*(E1sum.^2)/N);
sig2(sig2<0)=0;
Ebr = E1sum/N;
% $$$ Hm0 = 4*sqrt(sig2);
% $$$ mask= mask/N;
