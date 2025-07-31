function [visc,mask] = get_funwave_viscocity_breaking_1D(runNAME,rawDIR,dx,dy,Navg);
%
% USAGE: [visc,mask] = get_funwave_viscocity_breaking_1D(runNAME,rawDIR,dx,dy,Navg);
%
% runNAME: used for archiving stats
% rawDIR:  location of raw ascii data files
% (dy,dx): grid spacing (meters)
% visc:    hovmoller plot of alongshore maximum "nubrk"
% mask:    time-average of logical (h+eta)>=0.1;


% get structure with all elevation data
fprintf('working on: %s \n',runNAME);
eta_files   = dir([rawDIR,filesep,'eta_*']);
nubrk_files = dir([rawDIR,filesep,'nubrk_*']);
Nf          = length(nubrk_files);
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
visc  = zeros(Nf,Nx);
E2sum = zeros(1,Nx);
mask  = zeros(1,Nx);
%
% loop over Nf files, load and compute stats
for ii = 1:Nf
    % load the eta_ file
    fprintf('\tloading files associated with: %s \n',eta_files(ii).name);
    eta      = load([rawDIR,filesep,eta_files(ii).name],'-ascii');
    nubrk    = load([rawDIR,filesep,nubrk_files(ii).name],'-ascii');
    %
    % water depth
    H        = h+eta;
    H(H<=0)  = 0.001;
    %
    H_avg     = mean(H,1);
    mask      = mask  + double(H_avg>=0.1);
    visc(ii,:)= max(nubrk,[],1);
    N         = N     + 1;
end
mask= mask/N;

