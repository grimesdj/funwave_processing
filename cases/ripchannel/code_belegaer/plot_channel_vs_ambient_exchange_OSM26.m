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
SPREADs   = {'s00','s10'};
lineStyle = {'-',':'};
N = 3;

% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.5;
% $$$ ppos11 = [xm       ym         pw ph];
% $$$ ppos21 = [xm       ym+ph+ag   pw ph];
% $$$ ppos12 = [xm+pw+ag ym         pw ph];
% $$$ ppos22 = [xm+pw+ag ym+ph+ag   pw ph];
% $$$ ps    = [2*xm+2*(pw  2*ym+ag+2*ph];
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cbpos  = [xm+pw+ag ym       ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+2*(ag+ph)];
fig0 = figure('units','centimeters');
fig0.Position(3:4)=ps;
set(fig0,'papersize',ps,'paperposition',[0 0 ps]);
ax1 = axes('units','centimeters','position',ppos1);
ax2 = axes('units','centimeters','position',ppos2);
%
%
cm = colormap(cmocean('thermal',N+1));
cm = cm(1:N,:);
for jj=1:length(SPREADs)
    NAME = ['highRip-',BATHY,'-',SPREADs{jj}];
    dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat']);
    axes(ax2)
    hold on
    for kk = 1:N
        plot((dat.x-50)/(dat.bar_location(kk)-50), dat.Uex_mean_channel(:,kk), lineStyle{jj},'linewidth',2,'color',cm(kk,:))
    end

    axes(ax1)
    hold on
    for kk = 1:N
        plot((dat.x-50)/(dat.bar_location(kk)-50), dat.Uex_eddy_channel(:,kk), lineStyle{jj},'linewidth',2,'color',cm(kk,:))
    end
end
pLeg = plot(xlim,-999*[1 1],'-k',xlim,-999*[1 1],':k','linewidth',2)
legend(pLeg,{'$\sigma_\theta=0$','$\sigma_\theta=10$'},'interpreter','latex','fontsize',8)

ylabel(ax2,'$U_\mathrm{ex}$ [m/s]')
set(ax2,'xlim',[0 2.5],'ylim',[0 0.125],'tickdir','out','ticklabelinterpreter','latex','fontsize',10,'xticklabel',[])

ylabel(ax1,'$U_\mathrm{ex}$ [m/s]')
xlabel(ax1,'$(x-x_{sl})/x_b$ []')
set(ax1,'xlim',[0 2.5],'ylim',[0 0.125],'tickdir','out','ticklabelinterpreter','latex','fontsize',10)
grid([ax1 ax2],'on')

cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',num2str(dat.height'),'yaxislocation','right','xaxislocation','top','ydir','normal','xtick',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8,'ticklength',4*get(cb,'ticklength'))
xlabel(cb,'$H_0$ [m]','interpreter','latex','horizontalalignment','center')


figname = [figDIR,filesep,'Uex_vs_height_and_spread_',BATHY,'.pdf'];
exportgraphics(fig0,figname)
