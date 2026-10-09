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
BATHY     = 'barRip0';
lineStyle = {'-',':'};
N = 3;

% figure parameters
xm = 2.5;
ym = 2.5;
pw = 6;
ph = 4;
ag = 0.5;
% $$$ ppos11 = [xm       ym         pw ph];
% $$$ ppos21 = [xm       ym+ph+ag   pw ph];
% $$$ ppos12 = [xm+pw+ag ym         pw ph];
% $$$ ppos22 = [xm+pw+ag ym+ph+ag   pw ph];
% $$$ ps    = [2*xm+2*(pw  2*ym+ag+2*ph];
ppos1 = [xm       ym         pw ph];
% $$$ ppos2 = [xm       ym+ph+ag   pw ph];
cbpos  = [xm+pw+ag ym       ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+(ag+ph)];
fig0 = figure('units','centimeters');
fig0.Position(3:4)=ps;
set(fig0,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos1);
%
%
cm = colormap(cmocean('phase',N+2));
cm = cm(2:N+1,:);
NAME = ['spreadRip-',BATHY];
dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat']);


iX  = find(dat.x>=mean(dat.Xke_max),1,'first');

axes(ax1)
hold on,
plot(dat.spread,dat.Uex_channel(iX,:)     ,'o','markerfacecolor',cm(1,:),'markeredgecolor',cm(1,:)); hold on
plot(dat.spread,dat.Uex_mean_channel(iX,:),'s','markerfacecolor',cm(2,:),'markeredgecolor',cm(2,:)); hold on
plot(dat.spread,dat.Uex_eddy_channel(iX,:),'d','markerfacecolor',cm(3,:),'markeredgecolor',cm(3,:)); hold on

legend(ax1.Children([3,2,1]),{'$\langle{u}\rangle + \bar{u}$','$\langle{u}\rangle$','$\bar{u}$'},'interpreter','latex','location','southeast','fontsize',8)
grid([ax1],'on')

ylabel(ax1,'$U_\mathrm{ex}$ [m/s]')
xlabel(ax1,'$\sigma_\theta$ [$^\circ$]')
set(ax1,'xlim',[0 21],'ylim',[0 0.1],'tickdir','out','ticklabelinterpreter','latex','fontsize',10)

% $$$ cb = axes('units','centimeters','position',cbpos);
% $$$ imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(dat.height'),'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
% $$$ xlabel(cb,'$H_0$ [m]','interpreter','latex','horizontalalignment','center')


figname = [figDIR,filesep,'Uex_vs_spread_',BATHY,'.pdf'];
exportgraphics(fig0,figname)
