function [Fbx, Fbx_avg] = get_funwave_BreakingDissipation_1D(runNAME,rawDIR,dx,dy,Navg);
%
% USAGE: [Fbx, Fbx_avg] = get_funwave_BreakingDissipation_1D(runNAME,rawDIR,dx,dy,Navg);
%
% runNAME: used for archiving stats
% rawDIR:  location of raw ascii data files
% (dy,dx): grid spacing (meters)
% (Fbx, Fby):  viscosity breaking force
% u_avg:   backward average over Navg points
% Fbr:     1/(h+eta)* div( nubrk*(h+eta)*grad( u ) )


% get structure with all elevation data
fprintf('working on BrkDis for: %s \n',runNAME);
fprintf('!!!! Not loading BrkSrcY !!!!!!!\n')
% $$$ eta_files   = dir([rawDIR,filesep,'eta_*']);
% $$$ nubrk_files = dir([rawDIR,filesep,'nubrk_*']);
% $$$ p_files     = dir([rawDIR,filesep,'p_*']);
% $$$ q_files     = dir([rawDIR,filesep,'q_*']);
Brk_avg_files   = dir([rawDIR,filesep,'BrkDissX_*']);
Nf_avg          = length(Brk_avg_files);
%
Brk_files   = dir([rawDIR,filesep,'BrkSrcX_*']);
Nf          = length(Brk_files);
%
% get bathymetry file for filtering by water-depth
h      = load([rawDIR,filesep,'dep.out'],'-ascii');
%h      = 0.25*( h(1:end-1,1:end-1) + h(2:end,1:end-1) +...
%		h(2:end  ,1:end-1) + h(2:end,2:end) );
h1yavg = mean(h,1);
[Ny,Nx]= size(h);
%
% pre-allocate variables
%
% loop over Nf files, load and compute stats
for ii = 1:min(Nf,Nf_avg)
    % load the eta_ file
    brksrc      = load([rawDIR,filesep,Brk_files(ii).name],'-ascii');
    brksrc_avg  = load([rawDIR,filesep,Brk_avg_files(ii).name],'-ascii');    
    % 
    Fbx(ii,:)   = mean(brksrc,1);
    Fbx_avg(ii,:)   = mean(brksrc_avg,1);    
end
% $$$ Hm0 = 4*sqrt(sig2);
% $$$ mask= mask/N;
