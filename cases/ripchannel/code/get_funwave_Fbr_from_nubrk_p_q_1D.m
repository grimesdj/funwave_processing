function [Fbx, Fby] = get_funwave_Fbr_from_nubrk_p_q_1D(runNAME,rawDIR,dx,dy,Navg);
%
% USAGE: [Fbx, Fby] = get_funwave_Fbr_from_nubrk_p_q(runNAME,rawDIR,dx,dy,Navg);
%
% runNAME: used for archiving stats
% rawDIR:  location of raw ascii data files
% (dy,dx): grid spacing (meters)
% (Fbx, Fby):  viscosity breaking force
% u_avg:   backward average over Navg points
% Fbr:     1/(h+eta)* div( nubrk*(h+eta)*grad( u ) )


% get structure with all elevation data
fprintf('working on Fbr for: %s \n',runNAME);
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
    Fbx(ii,:)             = mean((Pyy+Pxx)./H,1); clear Pxx Pyy
    Fby(ii,:)             = mean((Qxx+Qyy)./H,1); clear Qxx Qyy
end
% $$$ Hm0 = 4*sqrt(sig2);
% $$$ mask= mask/N;
