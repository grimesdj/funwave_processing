function [Hs_xy,eta_bar,mask0,freqs,Snn_xy,xsl] = estimate_wave_stats_from_eta_h_t_dt(eta,h,t,dt,x);
%
% USAGE: [Hs_xy,eta_bar,mask0,freqs,Snn_xy,xsl] = estimate_wave_stats_from_eta_h_t_dt(eta,h,t,dt,x);
%

[ny,nx] = size(h);
nt = length(t);
% make a mask for eta
mask0  = (h + eta)>0.1;
% what should be the "wetted" threshold? (always?, 50%?, other?)
avgmask= mean(mask0,3);
mask0  = avgmask>0.01;% try 1% deeper than 0.1m!
% Make a shoreline mask were for all times h+eta>0.01
minmask = (h + min(eta,[],3))>0.01;
%
% first establish time-averaged waterlevel
eta_bar = nanmean(eta,3);
% convert eta from instantaneous waterlevel
% to perturbation waterlevel
eta     = (eta-eta_bar).*mask0;
eta_bar = eta_bar.*mask0;
%
% h_bar  = nanmean(h,1);
sig_eta=4*nanstd(eta,[],3);
%
% estimate the shoreline location (10cm depth; west coast/east coast)
dum = ones(ny,1)*[1:nx];
if h(1,1)<h(1,end)
% $$$     dum((h+eta_bar).*mask0-0.1<0)=inf;
    dum(minmask==0)=inf;
    ixsl = min(dum,[],2);
else
% $$$     dum((h+eta_bar).*mask0-0.1<0)=-inf;
    dum(minmask==0)=-inf;
    ixsl = max(dum,[],2);
end
xsl = x(ixsl);
%
% first reshape to have (t,[x,y])
dum = permute(eta,[3 1 2]);
dum2= reshape(dum, [nt, nx*ny]);
n2 = floor(log2(nt));
%
df0   = max(0.0125,1/(nt*dt));
chnk  = max(floor(nt*mean(dt)*df0),1);
inds  = 1:chnk*floor(1/(dt*df0));
eta_t_xy = dum2(inds,:);
[Snn,freqs] = welch_method(eta_t_xy,dt,chnk,0.0);
% some nans due to spetra of masked regions
Snn = real(Snn);
df      = freqs(2)-freqs(1);
fnquist = 0.5/dt;
nf      = size(Snn,1);
%
dum = reshape(Snn,nf,ny,nx);
dum2= nanmean(dum,2);
Snn_xy = dum;
Snn_x  = squeeze(dum2);
% estimate Hs
dum = sqrt(sum(16*Snn*df,1));
Hs_xy = reshape(dum,[ny,nx]);
Hs_x  = nanmean(Hs_xy,1);
%
mask0 = double(minmask);


