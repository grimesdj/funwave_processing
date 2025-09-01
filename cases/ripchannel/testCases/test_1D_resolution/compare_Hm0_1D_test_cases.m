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
ph = 3;
ag = 0.5;
ppos1 = [xm ym pw ph];
ps    = [2*xm+pw 2*ym+ph];
%
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
%
clrs = {'c','b','r'};
a1 = axes('units','centimeters','position',ppos1);
for ii=1:length(runBATHYlist)
    runBATHY = runBATHYlist{ii};
    runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
    matDIR   = [runDIR,filesep,'mat_data'];
    % load the original versions of Hm0
    load([matDIR,filesep,runBATHY,'_Hm0.mat'],'Hm0','mask','run_dirs','Ebr')
    Hm0_0 = (Hm0(end,:));
    Ebr_0 = Ebr(end,:);
    %
    % load the new versions of Hm0
    load([matDIR,filesep,runBATHY,'_test_case_Hm0.mat'],'Hm0','mask','run_dirs','Ebr')
    Hm0 = (Hm0(end,:));
    %
    N = size(Hm0,2);
    x = [0:N-1]'*dx(ii);
    if ii==1
        x_ref   = x;
        Hm0_ref = Hm0_0;
    end
    tmp = interp1(x_ref,Hm0_ref',x);
    tmp = tmp';
% $$$     whos tmp Hm0
    %
    axes(a1)
    hold(a1,'on')
    if ii>1
        plot(x,Hm0_0-tmp,'-','color',clrs{ii},'linewidth',1);
    end
    plt(ii) = plot(x,Hm0-tmp,'--','color',clrs{ii},'linewidth',2);
end
ylabel(a1,'$(1.5m,10s)$','interpreter','latex','fontsize',8)
xlabel(a1,'$x$ m','interpreter','latex','fontsize',8)
set(a1,'tickdir','out','ticklabelinterpreter','latex')%,'ylim',[0 1.75])
leg = legend(plt,'$dx=0.25$','$dx=0.50$m','$dx=1.00$m');
set(leg,'interpreter','latex','fontsize',7)
figname   = [figDIR,filesep,'Hm0_varying_resolution_test_cases.pdf'];
exportgraphics(fig,figname)
