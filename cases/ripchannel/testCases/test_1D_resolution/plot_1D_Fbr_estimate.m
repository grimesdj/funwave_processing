%% make a detailed plot of each term in Fbr.
% I'm looking to understand why Fbr switches sign
% mid-surfzone.
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))

% 1) determine which run to process
%
runBATHY = 'resolution1D';
figDIR   = '/storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/test_1D_resolution/figures/';
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
matDIR   = [runDIR,filesep,'mat_data'];
% $$$ % the list of run directories are saved in:
% $$$ % load([matDIR,filesep,'runs_to_process.mat'])
runNAME  = 'planar1Ddx050_h10t10s00d00';
rawDIR   = [runDIR,filesep,runNAME,filesep,'output'];

% 2) get the source for (eta,u,v,nubrk)
infoFile = [matDIR,filesep,'ripchannel_run_info_',runNAME,'.mat'];
[info,t0,dt0] = ripchannel_run_info_cmshpc(matDIR,runNAME);

BrkSrcFile = [matDIR,filesep,info.rootName,'BrkSrcX_01.nc'];
EtaFile    = [matDIR,filesep,info.rootName,'eta_01.nc'];
UFile      = [matDIR,filesep,info.rootName,'u_01.nc'];
NuBrkFile  = [matDIR,filesep,info.rootName,'nubrk_01.nc'];
DepFile    = [matDIR,filesep,info.rootName,'dep.nc'];

x = ncread(DepFile,'x');
dx= x(2)-x(1);
t = ncread(DepFile,'t');
h = ncread(DepFile,'dep');
eta     = ncread(EtaFile,'eta');
% create min depth mask
H    = h+eta;
mask = min(H+eta)>0.1;
%
BrkSrcX = ncread(BrkSrcFile,'BrkSrcX');
BrkSrcX(~mask)=nan;
% use rms to find time with strong breaking
rmsBrk  = rms(BrkSrcX,2,'omitnan');
[~,it]  = max(rmsBrk);
%
eta     = squeeze(eta(1,:,it));
H       = squeeze(H(1,:,it));
BrkSrcX = squeeze(BrkSrcX(1,:,it));
u       = ncread(UFile,'u',[1 1 it],[inf inf 1]);
nubrk   = ncread(NuBrkFile,'nubrk',[1 1 it],[inf inf 1]);
%
%% In funwave, nubrk should be equivalent to:
nubrk0  = H.*0.35.*sqrt(9.81*H);
nubrk0(nubrk==0)=0;% this is exactly correct!
%
% 3) estimate terms:
% dUHdx_fwd = (HU    (I+1,J)-HU    (I  ,J))/dx
% dUHdx_bwd = (HU    (I  ,J)-HU    (I-1,J))/dx
% nubrk_fwd = (nu_vis(I+1,J)+nu_vis(I  ,J))/2;
% nubrk_bwd = (nu_vis(I  ,J)+nu_vis(I-1,J))/2;
%
% Fbr = (nubrk_fwd*dUHdx_fwd - nubrk_bwd*dUHdx_bwd)/dx;
Fbr       = 0*H;
%
HU = u.*H;
dUHdx_cen = (HU(2:end)-HU(1:end-1))/dx;
nubrk_cen = (nubrk(1:end-1)+nubrk(2:end))/2;
Fbr(2:end-1) = (nubrk_cen(2:end).*dUHdx_cen(2:end)-nubrk_cen(1:end-1).*dUHdx_cen(1:end-1))/dx;
% 4) make 4-panel plot:
%   a) eta(x), 
%   b) uH(x), 
%   c) nubrk
%   d) Fbr, BrkSrcX
xlims = [100 125];
xm = 2.5;
ym = 2.5;
ag = 0.4;
pw = 10;
ph = 3.5;
ppos1 = [xm ym           pw ph];
ppos2 = [xm ym+ph+ag     pw ph];
ppos3 = [xm ym+2*(ph+ag) pw ph];
ppos4 = [xm ym+3*(ph+ag) pw ph];
ps    = [2*xm+pw 2*ym+4*ph+3*ag];
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize = ps;
fig.PaperPosition = [0 0 ps];
%
ax1 = axes('units','centimeters','position',ppos4);
plot(x,eta,'.k','linewidth',2)
grid on
ylabel(ax1,'$\eta$ [m]','interpreter','latex')
set(ax1,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax2 = axes('units','centimeters','position',ppos3);
plot(x,HU,'.k','linewidth',2)
grid on
ylabel(ax2,'$uH$ [m$^2$/s]','interpreter','latex')
set(ax2,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax3 = axes('units','centimeters','position',ppos2);
plot(x,nubrk,'.k','linewidth',2)
grid on
ylabel(ax3,'$\nu_\mathrm{br}$ [m$^2$/s]','interpreter','latex')
set(ax3,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax4 = axes('units','centimeters','position',ppos1);
plot(x,Fbr,'.k',x,BrkSrcX,'or','linewidth',2)
grid on
legend({'$F_{\nu}$','$F_\mathrm{tvd}$'},'interpreter','latex')
ylabel(ax4,'$F_\mathrm{br}$ (m/s)$^2$','interpreter','latex')
xlabel(ax4,'$x$ [m]','interpreter','latex')
set(ax4,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex')

figname = [figDIR,'breaking_force_two_ways_zoomed.pdf'];
exportgraphics(fig,figname)
%
%
%
%
%
%
eta     = ncread(EtaFile,'eta');
% create min depth mask
H    = h+eta;
mask = min(H+eta,[],3)>0.1;
%
BrkSrcX = ncread(BrkSrcFile,'BrkSrcX');
%
u       = ncread(UFile,'u');
nubrk   = ncread(NuBrkFile,'nubrk');
%
%
% 3) estimate terms:
% dUHdx_fwd = (HU    (I+1,J)-HU    (I  ,J))/dx
% dUHdx_bwd = (HU    (I  ,J)-HU    (I-1,J))/dx
% nubrk_fwd = (nu_vis(I+1,J)+nu_vis(I  ,J))/2;
% nubrk_bwd = (nu_vis(I  ,J)+nu_vis(I-1,J))/2;
%
% Fbr = (nubrk_fwd*dUHdx_fwd - nubrk_bwd*dUHdx_bwd)/dx;
Fbr       = 0*eta;
%
HU = u.*H;
dUHdx_cen = (HU(:,2:end,:)-HU(:,1:end-1,:))/dx;
nubrk_cen = (nubrk(:,1:end-1,:)+nubrk(:,2:end,:))/2;
Fbr(1,2:end-1,:) = (nubrk_cen(:,2:end,:).*dUHdx_cen(:,2:end,:)-nubrk_cen(:,1:end-1,:).*dUHdx_cen(:,1:end-1,:))/dx;
%
%
% now load the momentum BrkDissX
BrkDissFile = [matDIR,filesep,info.rootName,'MomentumTerms_01.nc'];
BrkDissX    = ncread(BrkDissFile,'BrkDissX');
%
% time average:
avgFbr    = mean(Fbr,3,'omitnan'); 
avgBrkSrc = mean(BrkSrcX,3,'omitnan');
avgBrkDiss= mean(BrkDissX,3,'omitnan');
%
% spatial average:
Nf = round(25/dx); if ~mod(Nf,2),Nf=Nf+1; end
flt = hamming(Nf); flt = flt./sum(flt);
avgFbr    = conv(avgFbr,flt,'same');
avgBrkSrc = conv(avgBrkSrc,flt,'same');
avgBrkDiss= conv(avgBrkDiss,flt,'same');
%
%
xlims = [50 200];
ppos  = [xm ym pw 2*ph];
ps    = [2*xm+pw 2*ym+2*ph];
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize = ps;
fig.PaperPosition = [0 0 ps];
%
ax1 = axes('units','centimeters','position',ppos);
plot(x,avgFbr,'-k',x,avgBrkSrc,'--r',x,avgBrkDiss,':b','linewidth',2)
%
legend({'$\bar{F}_{\nu}$','$\bar{F}_\mathrm{tvd}$','$\langle{F}_\mathrm{tvd}\rangle$'},'interpreter','latex')
ylabel(ax1,'$F_\mathrm{br}$ (m/s)$^2$','interpreter','latex')
xlabel(ax1,'$x$ [m]','interpreter','latex')
set(ax1,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex')
%
figname = [figDIR,'breaking_force_time_averaged.pdf'];
exportgraphics(fig,figname)
