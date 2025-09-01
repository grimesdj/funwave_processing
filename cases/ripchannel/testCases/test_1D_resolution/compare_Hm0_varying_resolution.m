% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
%
runBATHYlist = {'barred1Ddx025','barred1Ddx050','barred1DSWE'};%'planar1D'
dx   = [0.25,0.5,1];
dy   = 1;
%
%
xm = 2;
ym = 2;
pw = 8;
ph = 2;
ag = 0.5;
ppos6 = [xm ym pw ph];
ppos5 = [xm ym+ph+ag pw ph];
ppos4 = [xm ym+2*(ph+ag) pw ph];
ppos3 = [xm ym+3*(ph+ag) pw ph];
ppos2 = [xm ym+4*(ph+ag) pw ph];
ppos1 = [xm ym+5*(ph+ag) pw ph];
ps    = [2*xm+pw 2*ym+6*ph+5*ag];
%
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
%
clrs = {'b','r'};
a1 = axes('units','centimeters','position',ppos1);
a2 = axes('units','centimeters','position',ppos2);
a3 = axes('units','centimeters','position',ppos3);
a4 = axes('units','centimeters','position',ppos4);
a5 = axes('units','centimeters','position',ppos5);
a6 = axes('units','centimeters','position',ppos6);
for ii=1:length(runBATHYlist)
    runBATHY = runBATHYlist{ii};
    runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
    matDIR   = [runDIR,filesep,'mat_data'];
    load([matDIR,filesep,runBATHY,'_Hm0.mat'],'Hm0','mask','run_dirs','Ebr')
    N = size(Hm0,2);
    x = [0:N-1]*dx(ii);
    if ii==1
        x_ref   = x;
        Hm0_ref = Hm0;
        continue
    end
    tmp = interp1(x_ref,Hm0_ref',x);
    tmp = tmp';
% $$$     whos tmp Hm0
    %
    axes(a1)
    hold(a1,'on')
    plot(x,Hm0(1,:)-tmp(1,:),'color',clrs{ii-1},'linewidth',2);
    axes(a2)
    hold(a2,'on')
    plot(x,Hm0(2,:)-tmp(2,:),'color',clrs{ii-1},'linewidth',2);
    axes(a3)
    hold(a3,'on')
    plot(x,Hm0(3,:)-tmp(3,:),'color',clrs{ii-1},'linewidth',2);
    axes(a4)
    hold(a4,'on')
    plot(x,Hm0(4,:)-tmp(4,:),'color',clrs{ii-1},'linewidth',2);
    axes(a5)
    hold(a5,'on')
    plot(x,Hm0(5,:)-tmp(5,:),'color',clrs{ii-1},'linewidth',2);
    axes(a6)
    hold(a6,'on')
    plt(ii-1) = plot(x,Hm0(6,:)-tmp(6,:),'color',clrs{ii-1},'linewidth',2);
end
ylabel(a1,'$(0.5m,8s)$','interpreter','latex','fontsize',8)
set(a1,'tickdir','out','xticklabel',[],'ticklabelinterpreter','latex')%,'ylim',[0 0.75])
ylabel(a2,'$(0.5m,10s)$','interpreter','latex','fontsize',8)
set(a2,'tickdir','out','xticklabel',[],'ticklabelinterpreter','latex')%,'ylim',[0 0.75])
ylabel(a3,'$(1.0m,8s)$','interpreter','latex','fontsize',8)
set(a3,'tickdir','out','xticklabel',[],'ticklabelinterpreter','latex')%,'ylim',[0 1.25])
ylabel(a4,'$(1.0m,10s)$','interpreter','latex','fontsize',8)
set(a4,'tickdir','out','xticklabel',[],'ticklabelinterpreter','latex')%,'ylim',[0 1.25])
ylabel(a5,'$(1.5m,8s)$','interpreter','latex','fontsize',8)
set(a5,'tickdir','out','xticklabel',[],'ticklabelinterpreter','latex')%,'ylim',[0 1.75])
ylabel(a6,'$(1.5m,10s)$','interpreter','latex','fontsize',8)
set(a6,'tickdir','out','ticklabelinterpreter','latex')%,'ylim',[0 1.75])
xlabel(a6,'$x$ m','interpreter','latex','fontsize',8)
leg = legend(plt,'$dx=0.50$m','$dx=1.00$m');
set(leg,'interpreter','latex','fontsize',7)
figname   = '/scratch/grimesdj/ripchannel/Hm0_varying_resolution.pdf';
exportgraphics(fig,figname)
