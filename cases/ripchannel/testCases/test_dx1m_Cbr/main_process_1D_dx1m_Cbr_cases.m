% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% requires the input bathymetry name as top-dir
rootDIR = '/scratch/grimesdj/ripchannel/barred1Ddx1m_testCases/';
runDIRs = {'barred1D_Cbr050','barred1D_Cbr075','barred1D_Cbr100','barred1D_Cbr125','barred1D_Cbr150'};%'planar1D'
matDIR  = [rootDIR,filesep,'mat_data/'];
figDIR  = '/storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/test_dx1m_Cbr/figures/';
%
dx   = 1;
dy   = 1;
reprocess = 1;
%
% make hovmoller plots of breaking
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
clear Hm0 mask visc Ebr viscBrk
Navg=20;
if reprocess
for ii=1:length(runDIRs)
    %
    runNAME = runDIRs{ii}
    runDIR   = [rootDIR,runNAME];
    %
    %
    rawDIR = [runDIR,filesep,'output/'];
    [Hm0_tmp, mask_tmp]   = get_funwave_Hm0_from_eta(runNAME,rawDIR);
    [visc,~]      = get_funwave_viscosity_breaking_1D(runNAME,rawDIR);
    [viscBrk ,~]  = get_funwave_Fbr_from_nubrk_p_q_1D(runNAME,rawDIR,dx,dy,Navg);
    [Brk ,BrkAvg] = get_funwave_BreakingDissipation_1D(runNAME,rawDIR,dx,dy,Navg);    
    Ebr_tmp       = get_funwave_Ebr_from_nubrk_p_q(runNAME,rawDIR,dx,dy,Navg);
    Hm0(ii,:)     = Hm0_tmp;
    mask(ii,:)    = mask_tmp;
    Ebr(ii,:)     = Ebr_tmp;
    clf(fig0)
    a1 = axes('units','centimeters','position',ppos);
    % assume x = 1m resolution starting at x=0
    x = [0:(size(Hm0,2)-1)]*dx;
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
    imagesc(t,x,log10(abs(viscBrk)'))
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
    % create alongshore average and rms
    Fbr_nubrk_avg(ii,:) = mean(viscBrk,1   ,'omitnan');
    Fbr_nubrk_std(ii,:) = std (viscBrk,[],1,'omitnan');    
    %
    Fbr_BrkSrc_avg(ii,:) = mean(Brk,1   ,'omitnan');
    Fbr_BrkSrc_std(ii,:) = std (Brk,[],1,'omitnan');
    %
    Fbr_BrkAvg_avg(ii,:) = mean(BrkAvg,1   ,'omitnan');
    Fbr_BrkAvg_std(ii,:) = std (BrkAvg,[],1,'omitnan');    
    %    
    save([matDIR,filesep,runNAME,'_test_case_Hm0.mat'],'Hm0_tmp','mask_tmp','runDIRs','Ebr_tmp','viscBrk','Brk','BrkAvg')
end
    save([matDIR,filesep,'Cbrk_test_cases_Hm0.mat'],'Hm0','mask','runDIRs','Ebr','Fbr_nubrk_avg','Fbr_nubrk_std','Fbr_BrkSrc_avg','Fbr_BrkSrc_std','Fbr_BrkAvg_avg','Fbr_BrkAvg_std')
    %
else
    load([matDIR,filesep,'Cbrk_test_cases_Hm0.mat'],'Hm0','mask','runDIRs','Ebr','Fbr_nubrk_avg','Fbr_nubrk_std','Fbr_BrkSrc_avg','Fbr_BrkSrc_std','Fbr_BrkAvg_avg','Fbr_BrkAvg_std')
end
%
% make a figure, too.
xm = 2.5;
ym = 2;
pw = 8;
ph = 4;
ppos = [xm ym pw ph];
lpos = [xm+pw+0.05 ym 3 ph];
ps   = [2*xm+pw+3 1.5*ym+ph];
%
Ndirs = length(runDIRs);
clrs = prism(Ndirs);
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
colororder(fig, clrs);
ax  = axes('units','centimeters','position',ppos);
colororder(ax,clrs);
%
% assume x = 1m resolution starting at x=0
x = [0:(size(Hm0,2)-1)]*dx;
plot(x,Hm0,'-','linewidth',1.5)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$H_{m_0}$ [m]','interpreter','latex')
set(ax,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
ax2  = axes('units','centimeters','position',lpos);
set(ax2,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
figname = [figDIR,filesep,'test_case_Hm0.pdf'];
exportgraphics(fig,figname)
%
clf(fig)
colororder(fig, clrs);
ax  = axes('units','centimeters','position',ppos);
colororder(ax,clrs);
%
Nflt= round(5/dx); if ~mod(Nflt,2), Nflt=Nflt+1; end
flt = hamming(Nflt); flt = flt./sum(flt);
Ebr_lp = conv2(Ebr',flt,'same');
plot(x,Ebr_lp,'-','linewidth',1.5)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\epsilon_\mathrm{br}$ [m$^2$/s$^{3}$]','interpreter','latex')
set(ax,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
ax1  = axes('units','centimeters','position',lpos);
set(ax1,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
figname = [figDIR,filesep,'test_case_Ebr.pdf'];
exportgraphics(fig,figname)
%
%
%
%
clf(fig)
colororder(fig, clrs);
%
ax  = axes('units','centimeters','position',ppos);
colororder(ax,clrs);
plot(x,Fbr_nubrk_avg,'-','linewidth',1.5)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\epsilon_\mathrm{br}$ [m$^2$/s$^{3}$]','interpreter','latex')
set(ax,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
ax1  = axes('units','centimeters','position',lpos);
set(ax1,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
figname = [figDIR,filesep,'test_case_Ebr.pdf'];
exportgraphics(fig,figname)
close all
%
%
ppos = [xm ym pw ph];
ppos1= [xm ym+ph+0.5 pw ph];
lpos = [xm+pw+0.05 ym 3 ph];
ps   = [2*xm+pw+3 1.5*ym+2*ph+0.5];
%
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps]);
colororder(fig, clrs);
%
ax0  = axes('units','centimeters','position',ppos);
colororder(ax0,clrs);
plot(ax0,x,Fbr_nubrk_avg,'-','linewidth',1.5)
xlabel(ax0,'$x$ [m]','interpreter','latex')
ylabel(ax0,'$\bar{F}_\mathrm{br}$ [m/s$^{2}$]','interpreter','latex')
set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
axl  = axes('units','centimeters','position',lpos);
set(axl,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
%
%
ax1  = axes('units','centimeters','position',ppos1);
colororder(ax1,clrs);
plot(ax1,x,Fbr_nubrk_std,'-','linewidth',1.5)
xlabel(ax1,'$x$ [m]','interpreter','latex')
ylabel(ax1,'$\mathrm{std}(F_\mathrm{br})$ [m/s$^{2}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
figname = [figDIR,filesep,'test_case_Fbr_nubrk.pdf'];
exportgraphics(fig,figname)
%
%
%
clf(fig)
ax0  = axes('units','centimeters','position',ppos);
colororder(ax0,clrs);
plot(ax0,x,Fbr_BrkSrc_avg,'-','linewidth',1.5)
xlabel(ax0,'$x$ [m]','interpreter','latex')
ylabel(ax0,'$\bar{F}_\mathrm{br}$ [m/s$^{2}$]','interpreter','latex')
set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
axl  = axes('units','centimeters','position',lpos);
set(axl,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
%
%
ax1 = axes('units','centimeters','position',ppos1);
colororder(ax1,clrs);
plot(ax1,x,Fbr_BrkSrc_std,'-','linewidth',1.5)
xlabel(ax1,'$x$ [m]','interpreter','latex')
ylabel(ax1,'$\mathrm{std}(F_\mathrm{br})$ [m/s$^{2}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
figname = [figDIR,filesep,'test_case_Fbr_BrkSrcX.pdf'];
exportgraphics(fig,figname)
%
%
clf(fig)
ax0 = axes('units','centimeters','position',ppos);
colororder(ax0,clrs);
plot(ax0,x,Fbr_BrkAvg_avg,'-','linewidth',1.5)
xlabel(ax0,'$x$ [m]','interpreter','latex')
ylabel(ax0,'$\bar{F}_\mathrm{br}$ [m/s$^{2}$]','interpreter','latex')
set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
%
axl  = axes('units','centimeters','position',lpos);
set(axl,'ylim',[0 Ndirs+1],'xlim',[0 1],'xtick',[],'ytick',[],'ycolor','none','xcolor','none')
for ii=1:Ndirs
    % parse the wave height/perios/dir/spread info
    Cbr = regexp(runDIRs{ii},'(?<=Cbr)(..)','match');
    str    = sprintf('$C_\\mathrm{br}\\times%1.2f$',str2num(Cbr{1})/10);
    text(0,Ndirs+1-ii,str,'color',clrs(ii,:),'interpreter','latex')
end
%
%
ax1  = axes('units','centimeters','position',ppos1);
colororder(ax1,clrs);
plot(ax1,x,Fbr_BrkAvg_std,'-','linewidth',1.5)
xlabel(ax1,'$x$ [m]','interpreter','latex')
ylabel(ax1,'$\mathrm{std}(F_\mathrm{br})$ [m/s$^{2}$]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex','fontsize',15,'xlim',[0 500])
figname = [figDIR,filesep,'test_case_Fbr_BrkAvgX.pdf'];
exportgraphics(fig,figname)
%
close all


