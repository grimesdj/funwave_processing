%% make a detailed plot of each term in Fbr.
% I'm looking to understand why Fbr switches sign within the wave face
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
t = t(1:size(eta,3));

% create min depth mask
H    = h+eta;
mask = min(H+eta)>0.1;
%
BrkSrcX = ncread(BrkSrcFile,'BrkSrcX');
BrkSrcX(~mask)=nan;
% 
% 0) at two locations x0 = [75 125]
x0 = [75 125];
[~,ix] = min(abs(x-x0),[],1);
% 1) use zero-up crossing to identify ~100 waves
eta_avg = mean(eta,3);
zup         = squeeze((eta(1,ix,1:end-1)-eta_avg(1,ix)).*(eta(1,ix,2:end)-eta_avg(1,ix)));
zup(zup>0)  = 0;
zup(zup<0)  = 1;
zup( (eta(1,ix,1:end-1)-eta_avg(1,ix)) > 0 ) = 0;
zup = logical(zup);
% 2) extract cross-shore structure of each wave,
%    need the approximate wave-length of each zup wave
dt = mean(diff(t));
tp = t(1:end-1);
it1 = find(zup(1,:)');
it2 = find(zup(2,:)');
T1  = diff(tp(it1));
T2  = diff(tp(it2));
% 3) make an ensemble averaged wave shape (scale x using wavenumber given zero-up period)
xp = -1:0.01:1;
k1 = wavenumber_FunwaveTVD(2*pi./T1,H(ix(1)));
k2 = wavenumber_FunwaveTVD(2*pi./T2,H(ix(2)));
l1 = 2*pi./k1;
l2 = 2*pi./k2;

w1  = nan(length(T1),length(xp));
Fbr1= nan(length(T1),length(xp));
Nu1= nan(length(T1),length(xp));
HU1 = nan(length(T1),length(xp));
for jj=1:length(T1)
    it = it1(jj);%floor(0.5*(it1(jj)+it1(jj+1)));
    w1(jj,:) = interp1((x-x0(1))/l1(jj)',squeeze(eta(1,:,it)-eta_avg),xp);

    % load the current velocity and viscosity
    u       = ncread(UFile,'u',[1 1 it],[inf inf 1]);
    nubrk   = ncread(NuBrkFile,'nubrk',[1 1 it],[inf inf 1]);

    % create breaking force
    HH      = H(1,:,it);
    Fbr     = 0*HH;
    HU      = u.*HH;
  dUHdx_cen = (HU(2:end)-HU(1:end-1))/dx;
  nubrk_cen = (nubrk(1:end-1)+nubrk(2:end))/2;
  Fbr(2:end-1) = (nubrk_cen(2:end).*dUHdx_cen(2:end)-nubrk_cen(1:end-1).*dUHdx_cen(1:end-1))/dx;
  % archive
  Nu1(jj,:) = interp1((x-x0(1))/l1(jj)',nubrk,xp);  
  HU1(jj,:) = interp1((x-x0(1))/l1(jj)',HU,xp);
  Fbr1(jj,:) = interp1((x-x0(1))/l1(jj)',Fbr,xp);
end

w2  = nan(length(T2),length(xp));
Fbr2= nan(length(T2),length(xp));
Nu2= nan(length(T2),length(xp));
HU2 = nan(length(T2),length(xp));
for jj=1:length(T2)
    it = it2(jj);%floor(0.5*(it2(jj)+it2(jj+1)));
    w2(jj,:) = interp1((x-x0(2))/l2(jj)',squeeze(eta(1,:,it)-eta_avg),xp);

    % load the current velocity and viscosity
    u       = ncread(UFile,'u',[1 1 it],[inf inf 1]);
    nubrk   = ncread(NuBrkFile,'nubrk',[1 1 it],[inf inf 1]);

    % create breaking force
    HH      = H(1,:,it);
    Fbr     = 0*HH;
    HU      = u.*HH;
  dUHdx_cen = (HU(2:end)-HU(1:end-1))/dx;
  nubrk_cen = (nubrk(1:end-1)+nubrk(2:end))/2;
  Fbr(2:end-1) = (nubrk_cen(2:end).*dUHdx_cen(2:end)-nubrk_cen(1:end-1).*dUHdx_cen(1:end-1))/dx;
  % archive
  Nu2(jj,:) = interp1((x-x0(2))/l2(jj)',nubrk,xp);  
  HU2(jj,:) = interp1((x-x0(2))/l2(jj)',HU,xp);
  Fbr2(jj,:) = interp1((x-x0(2))/l2(jj)',Fbr,xp);
end

% 4) plot

% 4) make 8-panel plot:
%   a)/e) eta(x), 
%   b)/f) uH(x), 
%   c)/g) nubrk
%   d)/h) Fbr, BrkSrcX
xlims = [-1 1];
xm = 2.5;
ym = 2.5;
ag = 0.4;
pw = 8;
ph = 2.5;
ppos1 = [xm ym           pw ph];
ppos2 = [xm ym+ph+ag     pw ph];
ppos3 = [xm ym+2*(ph+ag) pw ph];
ppos4 = [xm ym+3*(ph+ag) pw ph];
% $$$ %
% $$$ ppos5 = [xm+pw+ag ym           pw ph];
% $$$ ppos6 = [xm+pw+ag ym+ph+ag     pw ph];
% $$$ ppos7 = [xm+pw+ag ym+2*(ph+ag) pw ph];
% $$$ ppos8 = [xm+pw+ag ym+3*(ph+ag) pw ph];
% $$$ %
% $$$ ps    = [2*xm+2*pw+ag 2*ym+4*ph+3*ag];
ps    = [2*xm+pw 2*ym+4*ph+3*ag];
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize = ps;
fig.PaperPosition = [0 0 ps];
%
ax1 = axes('units','centimeters','position',ppos4);
plot(xp,mean(w1,1),'-b',xp,mean(w2,1),'r','linewidth',2)
hold on,
plot(xp,mean(w1,1)+[-1; 1]*std(w1,0,1),':b',xp,mean(w2,1)+[-1; 1]*std(w2,0,1),':r','linewidth',2)
grid on
ylabel(ax1,'$\eta$ [m]','interpreter','latex')
set(ax1,'tickdir','out','fontsize',12,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax2 = axes('units','centimeters','position',ppos3);
plot(xp,mean(HU1,1),'-b',xp,mean(HU2,1),'-r','linewidth',2)
hold on,
plot(xp,mean(HU1,1)+[-1; 1]*std(HU1,0,1),':b',xp,mean(HU2,1)+[-1; 1]*std(HU2,0,1),':r','linewidth',2)
grid on
ylabel(ax2,'$uH$ [m$^2$/s]','interpreter','latex')
set(ax2,'tickdir','out','fontsize',12,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax3 = axes('units','centimeters','position',ppos2);
plot(xp,mean(Nu1,1,'omitnan'),'-b',xp,mean(Nu2,1,'omitnan'),'linewidth',2)
% $$$ Nu1(Nu1==0)=nan;
% $$$ Nu2(Nu2==0)=nan;
% $$$ plot(xp,10.^mean(log10(Nu1),1,'omitnan'),'-b',xp,10.^mean(log10(Nu2),1,'omitnan'),'linewidth',2)
% $$$ hold on,
% $$$ plot(xp,10.^(mean(log10(Nu1),1,'omitnan')+[-1; 1]*std(log10(Nu1),0,1,'omitnan')),':b',xp,10.^(mean(log10(Nu2),1,'omitnan')+[-1; 1]*std(log10(Nu2),0,1,'omitnan')),':r','linewidth',2)
grid on
ylabel(ax3,'$\nu_\mathrm{br}$ [m$^2$/s]','interpreter','latex')
set(ax3,'tickdir','out','fontsize',12,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[])
%
ax4 = axes('units','centimeters','position',ppos1);
plot(xp,mean(Fbr1,1),'-b',xp,mean(Fbr2,1),'-r','linewidth',2)
yline( mean( mean(Fbr1,1) ), '--b')
yline( mean( mean(Fbr2,1) ), '--r')
% $$$ hold on,
% $$$ plot(xp,mean(Fbr1,1)+[-1; 1]*std(Fbr1,0,1),':b',xp,mean(Fbr2,1)+[-1; 1]*std(Fbr2,0,1),'r','linewidth',2)
grid on
str = {['$x_0=',num2str(x0(1)),'$ m'],['$x_0=',num2str(x0(2)),'$ m']};
legend(ax1,str,'interpreter','latex','orientation','vertical','location','northwest','fontsize',6)
ylabel(ax4,'$F_\mathrm{br}$ (m/s)$^2$','interpreter','latex')
xlabel(ax4,'$(x-x_0)/\lambda$ [n/a]','interpreter','latex')
set(ax4,'tickdir','out','fontsize',12,'xlim',xlims,'ticklabelinterpreter','latex')
str = {['$\bar{F}_\mathrm{br} = ',num2str(mean(mean(Fbr1)),'%1.2e'),'$'], ['$\bar{F}_\mathrm{br} = ',num2str(mean(mean(Fbr2)),'%1.2e'),'$']};
legend(ax4,str,'interpreter','latex','location','southwest','fontsize',6)
% $$$ %% second column
% $$$ 
% $$$ ax5 = axes('units','centimeters','position',ppos8);
% $$$ plot(x,eta,'.k','linewidth',2)
% $$$ grid on
% $$$ set(ax5,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
% $$$ %
% $$$ ax6 = axes('units','centimeters','position',ppos7);
% $$$ plot(x,HU,'.k','linewidth',2)
% $$$ grid on
% $$$ %ylabel(ax2,'$uH$ [m$^2$/s]','interpreter','latex')
% $$$ set(ax6,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
% $$$ %
% $$$ ax7 = axes('units','centimeters','position',ppos6);
% $$$ plot(x,nubrk,'.k','linewidth',2)
% $$$ grid on
% $$$ %ylabel(ax3,'$\nu_\mathrm{br}$ [m$^2$/s]','interpreter','latex')
% $$$ set(ax7,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
% $$$ %
% $$$ ax8 = axes('units','centimeters','position',ppos5);
% $$$ plot(x,Fbr,'.k',x,BrkSrcX,'or','linewidth',2)
% $$$ grid on
% $$$ legend({'$F_{\nu}$','$F_\mathrm{tvd}$'},'interpreter','latex')
% $$$ ylabel(ax8,'$F_\mathrm{br}$ (m/s)$^2$','interpreter','latex')
% $$$ xlabel(ax8,'$x$ [m]','interpreter','latex')
% $$$ set(ax8,'tickdir','out','fontsize',15,'xlim',xlims,'ticklabelinterpreter','latex')

figname = [figDIR,'breaking_force_ensembles_inner_vs_outer_surfzone.pdf'];
exportgraphics(fig,figname)
