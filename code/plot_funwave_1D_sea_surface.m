function [eta_bar,sig_eta,eta_xt,h_bar] = plot_funwave_1D_sea_surface(rootMat,rootName,bathyFile,Hs,Tp,subDomain,plotter);
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
%
% load grid and bottom 
load(bathyFile);
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
files = dir([rootMat,rootName,'eta*.mat']);
eta0  = [];
t0    = [];
for ii=1:length(files)
    load([files(ii).folder,'/',files(ii).name]);
    if exist('subDomain','var')
        eta = eta(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
    end
    eta0 = cat(3,eta0,eta);
    t0   = cat(1,t0,t);
end
eta = eta0;
t   = t0;
dt  = mean(diff(t));
%
nt  = length(t);
h_bar = nanmean(h,1);
eta_bar=nanmean(nanmean(eta,3),1);
%
sig_eta=4*sqrt(mean(nanvar(eta,[],3),1));
%
