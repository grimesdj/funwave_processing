%% plot the eddy velocity at the center of the rip:
clear all
close all
%
%
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures/'
outDIR   = '/data2/ripchannel/mat_data/'
%
%
%% 1) need a list of run-directories
BATHYs    = {'barRip0','barRip1','terRip1'};
lineStyle = {'-','--',':'};
tau = [1000 700 725];
t0  = [1600 750 450];
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.5;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cbpos = [xm+pw+ag ym       ag ph/2];
ps    = [2*xm+pw  2*ym+ag+2*ph];
fig0 = figure('units','centimeters');
fig0.Position(3:4)=ps;
set(fig0,'papersize',ps,'paperposition',[0 0 ps]);
ax = axes('units','centimeters','position',ppos1);
for jj=1:length(BATHYs)
    NAME = ['spreadRip-',BATHYs{jj}];
    dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat']);
    hold on
% $$$     info = load([rootDIR,'/spreadRip/mat_data/ripchannel_run_info_',BATHYs{jj},'_h10t10s00d00.mat']);
% $$$     fin = [rootDIR,'/spreadRip/mat_data/funwave_',BATHYs{jj},'_h10t10s00d00_velocity_decomposition.nc'];
% $$$     Urot = ncread(fin,'Urot');
% $$$     Urot_mean = mean(Urot,3);
% $$$     Urot = Urot-Urot_mean;
% $$$     Urip = squeeze(Urot(round(info.Ly/2),round(info.xc/info.dx),:));
% $$$     plot(([1:30:30*200]-t0(jj))./tau(jj),Urip,lineStyle{jj},'color','k')
    plot(([1:30:30*200]-t0(jj))./tau(jj),dat.Umax_eddy_vs_t(:,1),lineStyle{jj},'color','k','linewidth',2)
    Urms(jj,:) = rms(dat.Umax_eddy_vs_t,1,'omitnan');
    Umean(jj,:)= dat.Umax_mean;
end
str = {sprintf('B(100-m): $\\tau=%d\\,\\mathrm{s}$',tau(1)),...
       sprintf('B(~50-m): $\\,\\tau=%d\\,\\mathrm{s}$',tau(2)),...
       sprintf('T(~50-m): $\\,\\tau=%d\\,\\mathrm{s}$',tau(3))};
legend(str,'interpreter','latex','fontsize',8,'units','centimeters','position',[xm-0.5+pw*2/3 ym+ph+ag pw/3 ph/2])
ylabel('$\bar{u}_\mathrm{rip}$ [m/s]')
xlabel('$(t-t_0)/\tau$ []')
grid on
set(ax,'xlim',[0 5],'tickdir','out','ticklabelinterpreter','latex','fontsize',10)

figname = [figDIR,filesep,'Urip_zero_spread_vs_bathy.pdf'];
exportgraphics(fig0,figname)
