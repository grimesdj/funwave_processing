function [Hm0,mask] = get_funwave_Hm0_from_eta(runNAME,rawDIR);
%
% USAGE: [Hm0,mask] = get_funwave_Hm0_from_eta(runNAME,rawDIR);
%
% runNAME: used for archiving stats
% rawDIR:  location of raw ascii data files
% Hm0:     4*std(eta) \approx 4*sqrt( mean(eta.^2)/(N-1) - mean(eta).^2)/(N*(N-1)) )
% mask:    fraction of time (h+eta)>0.1 m

% get structure with all elevation data
fprintf('working on Hm0 for: %s \n',runNAME);
files = dir([rawDIR,filesep,'eta_*']);
Nf    = length(files);
%
% get bathymetry file for filtering by water-depth
h     = load([rawDIR,filesep,'dep.out'],'-ascii');
h1yavg= mean(h,1);
[Ny,Nx]= size(h);
%
% pre-allocate variables
N     = 0;
e1sum = zeros(1,Nx);
e2sum = zeros(1,Nx);
mask  = zeros(1,Nx);
%
% loop over Nf files, load and compute stats
for ii = 1:Nf
    % load the eta_ file
    % fprintf('loading: %s \n',files(ii).name);    
    eta      = load([files(ii).folder,filesep,files(ii).name],'-ascii');
    eta1yavg = nanmean(eta,1);
    eta2yavg = nanmean(eta.^2,1);
    e1sum    = e1sum + eta1yavg;
    e2sum    = e2sum + eta2yavg;
    mask     = mask  + double(eta1yavg+h1yavg>=0.1);
    N        = N     + 1;
end
sig2= Ny/( (N-1)*(Ny-1)) * (e2sum - 2*(e1sum.^2)/N);
sig2(sig2<0)=0;
Hm0 = 4*sqrt(sig2);
mask= mask/N;
