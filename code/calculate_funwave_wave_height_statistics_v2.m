function [Hs_xy,eta_bar,mask0,x,y,h,freqs,Snn_xy,xsl] = calculate_funwave_wave_height_statistics_v2(rootMat,rootName,bathyFile,Hs,Tp,subDomain);
%
% USAGE: [Hs_x,x,eta_bar,h_bar,sig_eta,freqs,Snn_xy,Snn_x,Snn_wg,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(rootMat,rootName,bathyFile,Hs,Tp,subDomain,plotter,iwg);
%
% plot sea-surface elevation
% need to add some way of figuring out the shoreline location
% xsl = 498;% x=0 shoreline location
%
% this version also loads the land-mask and uses that to remove points near the "shoreline" to limit how much swash motions are analyze. 
%
% file_mat = [root,bath,'/',run,'/','dep.mat'];
% h = dep;

% load grid and bottom 
    load(bathyFile);
    h = -h;
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
files = dir([rootMat,rootName,'eta*.mat']);
% mask0  = [];
eta0   = [];
t0     = [];
for ii=1:length(files)
    fprintf('loading eta from: %s \n', files(ii).name);
    if exist('subDomain','var')
        dat0= matfile([files0(ii).folder,'/',files0(ii).name]);
        dat = matfile([files(ii).folder,'/',files(ii).name]);
        %        mask= dat0.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        eta = dat.eta(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        t = dat.t;
    else
        load([files0(ii).folder,'/',files0(ii).name]);        
        load([files(ii).folder,'/',files(ii).name]);
    end
    %    mask0 = cat(3,mask0,mask);    
    eta0 = cat(3,eta0,eta);
    t0   = cat(1,t0,t);
end
%mask = mask0;
eta = eta0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
%
clear eta0 t0 %mask0
%
% make a mask for eta
mask0  = (h + eta)>0.1;
% what should be the "wetted" threshold? (always?, 50%?, other?)
avgmask= mean(mask0,3);
mask0  = avgmask>0.01;% try 99% deeper than 0.1m!
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
% $$$ eta_bar=nanmean(eta,3);
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
%
% $$$ cm = cmocean('ballance');
% $$$ figure
% $$$ for ii=2:nt
% $$$     %    plot(x,squeeze(eta(128,:,ii)),'-k',x,squeeze(eta(128,:,ii)),':r')
% $$$     contourf(x,y,squeeze(eta(:,:,ii)),[-0.25:0.01:0.25],'edgecolor','none')
% $$$     caxis([-0.25 0.25]), colormap(cm)
% $$$     pause(0.1)
% $$$ end
%
% compute wave statistics
% first reshape to have (t,[x,y])
dum = permute(eta,[3 1 2]);
dum2= reshape(dum, [nt, nx*ny]);
n2 = floor(log2(nt));
% $$$ eta_t_xy = dum2(nt-2^n2+1:end,:);
% $$$ [Snn,nfft,df,fnquist,freqs,err_low,err_high] = Spectra(eta_t_xy,4,0.75,1,inf);
df0   = max(0.0125,1/(nt*dt));
chnk  = max(floor(nt*mean(dt)*df0),1);
inds  = 1:chnk*floor(1/(dt*df0));
eta_t_xy = dum2(inds,:);
[Snn,freqs] = mywelch(eta_t_xy,dt,chnk,0.0);
% some nans due to spetra of masked regions
Snn = real(Snn);
df      = freqs(2)-freqs(1);
fnquist = 0.5/dt;
nf      = size(Snn,1);
%
% $$$ % create synthetic Jonswap/PM spectrum
% $$$ if ~exist('Hs','var')
% $$$     Hs = 0.5;
% $$$     Tp = 8;
% $$$ end
% $$$ gam= 1;
% $$$ om = 2*pi*freqs;
% $$$ [S,~,~]=Jonswap('Omega', om ,'Hs',Hs,'Tm' ,Tp, 'Type', 2, 'TEnd',90,'Cap',3);
% $$$ SJONSWAP = 2*pi*S;
%
dum = reshape(Snn,nf,ny,nx);
dum2= nanmean(dum,2);
Snn_xy = dum;
Snn_x  = squeeze(dum2);
%
% $$$ % for funwave, look at wave-guage
% $$$ % statistics away from wave-maker
% $$$ if ~exist('iwg','var')
% $$$     %
% $$$     % need to be onshore of wave maker
% $$$     [~,iwg] = min(abs(h_bar-0.99*min(h_bar)));
% $$$ end
% $$$ Snn_wg = squeeze(Snn_xy(:,1,iwg));
%
%
% estimate Hs
dum = sqrt(sum(16*Snn*df,1));
Hs_xy = reshape(dum,[ny,nx]);
Hs_x  = nanmean(Hs_xy,1);
%
mask0 = minmask;
% $$$ % try using linear theory to shoal spectrum
% $$$ g = 9.8;
% $$$ k0 = wavenumber(om',-h_bar);
% $$$ dpth=ones(length(om),1)*max((eta_bar-h_bar),0.1);
% $$$ Cg = 0.5*(g*tanh(k0.*dpth)+g*k0.*dpth.*(sech(k0.*dpth).^2))./sqrt(g*k0.*tanh(k0.*dpth));
% $$$ Snn_linear = (Snn_wg.*Cg(:,iwg).*ones(1,length(h_bar)))./Cg;
% $$$ Hs_linear = sqrt(nansum(16*Snn_linear*df,1));
% $$$ %
% $$$ % this is for west-coast orientation
% $$$ is = min(iwg,ishoreline);
% $$$ ie = max(iwg,ishoreline);
% $$$ %
% $$$ % estimate location of break-point
% $$$ [Hsz,iHsz] = nanmax(Hs_x(is:ie));
% $$$ isurfzone  = iHsz+is;
% $$$ xsz = x(isurfzone);
% $$$ Lsz = xsl-xsz;
% $$$ %
% $$$ Snn_onshore = squeeze(dum2(:,1,isurfzone));
% $$$ %
% $$$ if plotter
% $$$     fig1=figure('color','w');
% $$$     I = find(freqs>0.05 & freqs<0.3);
% $$$     loglog(freqs(I), Snn_wg(I),'-k',freqs(I),S(I)*2*pi,'--r',freqs(I),Snn_onshore(I),'--b')
% $$$     set(gca,'ylim',[10^-4 10^0],'xlim',[1/128 1])
% $$$     ylims = ylim(gca);
% $$$     hold on,plot(1/Tp*[1 1],ylims)
% $$$     xlabel('$f$~Hz','interpreter','latex')
% $$$     ylabel('$S_{\eta\eta}~[\mathrm{m^2\,Hz^{-1}}]$','interpreter','latex')
% $$$ if ~exist([root,bath,'/','figures'],'dir')
% $$$     eval(['!mkdir ',root,bath,'/','figures/'])
% $$$ end
% $$$ set(gca,'fontsize',15)
% $$$ f1name = sprintf('%s%s/figures/Snn_wg_%s.pdf',root,bath,run);
% $$$ exportgraphics(fig1,f1name)
% $$$ %
% $$$ fig2 = figure('color','w'),
% $$$ ax = subplot(2,1,1);
% $$$ plot(x-xsl, Hs_x,'-k',xsz-xsl,Hsz,'xr')
% $$$ set(ax,'xlim',[-300 10],'ylim',[0 1.2*Hs])
% $$$ grid on
% $$$ % ax = get(gca);
% $$$ % $$$ xlabel('$x$','interpreter','latex')
% $$$ ylabel('$H_s(x)~\mathrm{[m]}$','interpreter','latex')
% $$$ ax2 = subplot(2,1,2);
% $$$ plot(x-xsl,-h_bar,'-k','linewidth',1)
% $$$ hold on, plot(x-xsl,eta_bar,'-b','linewidth',0.5)
% $$$ xlabel('$x$','interpreter','latex')
% $$$ ylabel('$h(x)~\mathrm{[m]}$')
% $$$ set(ax2,'xlim',ax.XAxis.Limits,'xtick',ax.XAxis.TickValues,'ylim',[-2.5 nanmax(eta_bar)])
% $$$ grid on
% $$$ % $$$ set(ax2,'color','none','Position',ax.Position,'InnerPosition',ax.InnerPosition,'yaxislocation','right',...
% $$$ % $$$         'xlim',ax.XAxis.Limits,'xtick',[])
% $$$ f2name = sprintf('%s%s/figures/Hs_vs_x_%s.pdf',root,bath,run);
% $$$ exportgraphics(fig2,f2name)
% $$$ %
% $$$ fig3=figure('color','w');
% $$$ cm = cmocean('balance');
% $$$ clrs = [0.75:0.5/255:1.25];
% $$$ contourf(x,y,sig_eta./repmat(Hs_x,ny,1)-1,clrs,'edgecolor','none')
% $$$ xlabel('$x$~[m]','interpreter','latex')
% $$$ ylabel('$y$~[m]','interpreter','latex')
% $$$ set(gca,'fontsize',15,'xlim',[-300 10],'color',0.5*[1 1 1])
% $$$ caxis([-0.25 0.25]),colormap(cm)
% $$$ cb = colorbar;
% $$$ ylabel(cb,'$H_s/\bar{H}_s-1$','interpreter','latex')
% $$$ f3name = sprintf('%s%s/figures/stdETAx4_%s.pdf',root,bath,run);
% $$$ exportgraphics(fig3,f3name)
% $$$ %
end