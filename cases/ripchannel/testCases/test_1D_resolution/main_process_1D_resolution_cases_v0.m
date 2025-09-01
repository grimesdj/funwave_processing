% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% requires the input bathymetry name as top-dir
runBATHYlist = {'barred1DSWE','barred1Ddx050','barred1Ddx025'};%'planar1D'
figDIR = '/storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/test_1D_resolution/figures/';
Navg = 20;% number of time-steps to average for breaking wave dissipation rate
dx   = [1,0.5,0.25];
dy   = 1;
for ii=1:length(runBATHYlist)
    runBATHY=runBATHYlist{ii};
    fprintf('working on: %s',runBATHY)
    %
    runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
    matDIR   = [runDIR,filesep,'mat_data'];
    %
    % the list of run directories are saved in:
    load([matDIR,filesep,'runs_to_process.mat'])
    % brings in cell array: run_dirs
    % for example,
    % run_dirs =
    %   6x1 cell array
    %    {'planar1D_h05t08s00d00'}
    %    {'planar1D_h05t10s00d00'}
    %    {'planar1D_h10t08s00d00'}
    %    {'planar1D_h10t10s00d00'}
    %    {'planar1D_h15t08s00d00'}
    %    {'planar1D_h15t10s00d00'}
    %
    % loop over run_dirs
    Ndirs  = length(run_dirs);
    %
    % create hovmoller plot of nubrk
    xm = 2.5;
    ym = 2;
    pw = 8;
    ph = 4;
    ppos = [xm ym pw ph];
    cpos = [xm+pw+0.05 ym 0.1 ph/2];
    ps   = [2*xm+pw+0.25 1.5*ym+ph];
    %
    fig0   = figure('units','centimeters');
    fig0.Position(3:4) = ps;
    fig0.PaperSize = ps;
    fig0.PaperPosition = [0 0 ps];
    %
    cm=cmocean('thermal');
    %
    clear Hm0 mask visc Ebr
    jj=Ndirs
    runNAME= run_dirs{jj};
    rawDIR = [runDIR,filesep,runNAME,filesep,'output'];
    [Hm0(jj,:), mask(jj,:)]   = get_funwave_Hm0_from_eta(runNAME,rawDIR);
    [visc,~]      = get_funwave_viscosity_breaking_1D(runNAME,rawDIR);
    [viscBrk ,~]  = get_funwave_Fbr_from_nubrk_p_q_1D(runNAME,rawDIR,dx(ii),dy,Navg);
    [Brk ,BrkAvg] = get_funwave_BreakingDissipation_1D(runNAME,rawDIR,dx(ii),dy,Navg);    
    Ebr(jj,:)                = get_funwave_Ebr_from_nubrk_p_q(runNAME,rawDIR,dx(ii),dy,Navg);
    clf(fig0)
    a1 = axes('units','centimeters','position',ppos);
    % assume x = 1m resolution starting at x=0
    x = [0:(size(Hm0,2)-1)]*dx(ii);
    t = 0:(size(visc,1)-1);
    imagesc(t,x,log10(visc'))
    clim = [-6 max(-5,log10(max(visc(:))))];
    clrs = [clim(1):diff(clim)/255:clim(2)];
    caxis(clim)
    colormap(cm)
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    title(runNAME)
    set(a1,'tickdir','out','ydir','normal','ticklabelinterpreter','latex')
    %
    c1 = axes('units','centimeters','position',cpos);
    imagesc(0,clrs,reshape(cm,256,1,3))
    set(c1,'yaxislocation','right','xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','xtick',[])
    xlabel('$\mathrm{log}_{10}\nu_\mathrm{br}$','interpreter','latex')
    figname = [figDIR,filesep,runNAME,'_test_case_nubrk.pdf'];
    exportgraphics(fig0,figname)
    %
    clf(fig0)
    a1 = axes('units','centimeters','position',ppos);
    % assume x = 1m resolution starting at x=0
    t = 0:(size(viscBrk,1)-1);
    imagesc(t,x,log10(abs(viscBrk')))
    clim = [-6 max(-5,log10(max(abs(viscBrk(:)))))];
    clrs = [clim(1):diff(clim)/255:clim(2)];
    caxis(clim)
    colormap(cm)
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    title(runNAME)
    set(a1,'tickdir','out','ydir','normal','ticklabelinterpreter','latex')
    %
    c1 = axes('units','centimeters','position',cpos);
    imagesc(0,clrs,reshape(cm,256,1,3))
    set(c1,'yaxislocation','right','xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','xtick',[])
    xlabel('$\mathrm{log}_{10}F_\mathrm{br}$','interpreter','latex')
    figname = [figDIR,filesep,runNAME,'_test_case_Fbr_from_nubrk.pdf'];
    exportgraphics(fig0,figname)
    %
    clf(fig0)
    a1 = axes('units','centimeters','position',ppos);
    % assume x = 1m resolution starting at x=0
    t = 0:(size(Brk,1)-1);
    imagesc(t,x,log10(abs(Brk')))
    clim = [-6 max(-5,log10(max(abs(Brk(:)))))];
    clrs = [clim(1):diff(clim)/255:clim(2)];
    caxis(clim)
    colormap(cm)
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    title(runNAME)
    set(a1,'tickdir','out','ydir','normal','ticklabelinterpreter','latex')
    %
    c1 = axes('units','centimeters','position',cpos);
    imagesc(0,clrs,reshape(cm,256,1,3))
    set(c1,'yaxislocation','right','xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','xtick',[])
    xlabel('$\mathrm{log}_{10}F_\mathrm{br}$','interpreter','latex')
    figname = [figDIR,filesep,runNAME,'_test_case_Fbr_from_BrkSrcX.pdf'];
    exportgraphics(fig0,figname)
    %
        clf(fig0)
    a1 = axes('units','centimeters','position',ppos);
    % assume x = 1m resolution starting at x=0
    t = 0:(size(BrkAvg,1)-1);
    imagesc(t,x,log10(abs(BrkAvg')))
    clim = [-6 max(-5,log10(max(abs(BrkAvg(:)))))];
    clrs = [clim(1):diff(clim)/255:clim(2)];
    caxis(clim)
    colormap(cm)
    xlabel('$t$ [s]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    title(runNAME)
    set(a1,'tickdir','out','ydir','normal','ticklabelinterpreter','latex')
    %
    c1 = axes('units','centimeters','position',cpos);
    imagesc(0,clrs,reshape(cm,256,1,3))
    set(c1,'yaxislocation','right','xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','xtick',[])
    xlabel('$\mathrm{log}_{10}F_\mathrm{br}$','interpreter','latex')
    figname = [figDIR,filesep,runNAME,'_test_case_Fbr_from_BrkDisX.pdf'];
    exportgraphics(fig0,figname)
    %
    save([matDIR,filesep,runBATHY,'_test_case_Hm0.mat'],'Hm0','mask','run_dirs','Ebr','viscBrk','Brk','BrkAvg')
% $$$ %
% $$$ %
% $$$ % make a figure, too.
% $$$ xm = 2.5;
% $$$ ym = 2;
% $$$ pw = 8;
% $$$ ph = 4;
% $$$ ppos = [xm ym pw ph];
% $$$ lpos = [xm+pw+0.05 ym 3 ph];
% $$$ ps   = [2*xm+pw+3 1.5*ym+ph];
% $$$ %
% $$$ clrs = prism(Ndirs);
% $$$ fig = figure('units','centimeters');
% $$$ pos = get(fig,'position');
% $$$ pos(3:4)=ps;
% $$$ set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
% $$$ colororder(fig, clrs);
% $$$ ax  = axes('units','centimeters','position',ppos);
% $$$ colororder(ax,clrs);
% $$$ %
% $$$ % assume x = 1m resolution starting at x=0
% $$$ x = [0:(size(Hm0,2)-1)]*dx(ii);
% $$$ plot(x,Hm0,'-','linewidth',1.5)
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$H_{m_0}$ [m]','interpreter','latex')
% $$$ set(ax,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
% $$$ %
% $$$ ax  = axes('units','centimeters','position',lpos);
% $$$ set(ax,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
% $$$ for jj=1:Ndirs
% $$$     % parse the wave height/perios/dir/spread info
% $$$     Hs = regexp(run_dirs{jj},'(?<=h)(..)','match');
% $$$     Tp = regexp(run_dirs{jj},'(?<=t)(..)','match');
% $$$     str    = sprintf('$H_\\mathrm{s}$=%1.1f, $T_\\mathrm{p}$=%d',str2num(Hs{1})/10,str2num(Tp{1}));
% $$$     text(0,Ndirs+1-jj,str,'color',clrs(jj,:),'interpreter','latex')
% $$$ end
% $$$ figname = [figDIR,filesep,runBATHY,'_test_case_Hm0.pdf'];
% $$$ exportgraphics(fig,figname)
%
%
% $$$ %
% $$$ %
% $$$ fig = figure('units','centimeters');
% $$$ pos = get(fig,'position');
% $$$ pos(3:4)=ps;
% $$$ set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
% $$$ colororder(fig, clrs);
% $$$ ax  = axes('units','centimeters','position',ppos);
% $$$ colororder(ax,clrs);
% $$$ %
% $$$ Nflt= round(5/dx(ii)); if ~mod(Nflt,2), Nflt=Nflt+1; end
% $$$ flt = hamming(Nflt); flt = flt./sum(flt);
% $$$ Ebr_lp = conv2(Ebr',flt,'same');
% $$$ plot(x,Ebr_lp,'-','linewidth',1.5)
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$\epsilon_\mathrm{br}$ [m$^2$/s$^{3}$]','interpreter','latex')
% $$$ set(ax,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
% $$$ %
% $$$ ax  = axes('units','centimeters','position',lpos);
% $$$ set(ax,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
% $$$ for jj=1:Ndirs
% $$$     % parse the wave height/perios/dir/spread info
% $$$     Hs = regexp(run_dirs{jj},'(?<=h)(..)','match');
% $$$     Tp = regexp(run_dirs{jj},'(?<=t)(..)','match');
% $$$     str    = sprintf('$H_\\mathrm{s}$=%1.1f, $T_\\mathrm{p}$=%d',str2num(Hs{1})/10,str2num(Tp{1}));
% $$$     text(0,Ndirs+1-jj,str,'color',clrs(jj,:),'interpreter','latex')
% $$$ end
% $$$ figname = [figDIR,filesep,runBATHY,'_test_case_Ebr.pdf'];
% $$$ exportgraphics(fig,figname)
% $$$ close all
end

compare_Hm0_1D_test_cases
