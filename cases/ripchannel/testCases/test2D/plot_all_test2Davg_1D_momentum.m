% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHY = 'test2Davg'
%
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
matDIR   = [runDIR,filesep,'mat_data'];
figDIR   = [runDIR,filesep,'figures/'];
%
% the list of run directories are saved in:
load([matDIR,filesep,'runs_to_process.mat'])
% run_dirs = cat(1,{'plnr2D_h10t10s10d00'},run_dirs);
%
% loop over run_dirs
Ndirs  = length(run_dirs);
for jj = 1:Ndirs
% 1) get current run subdirectory to process:
runID    = run_dirs{jj};
fprintf('\n processing: %s %s \n', runBATHY,runID)    
% 2) get the archived info structure:
if ismember(runID,'plnr2D_h10t10s10d00','rows')
    infoFile = dir([runDIR,filesep,'..',filesep,'test2D',filesep,'mat_data',filesep,'*','info','*',runID,'.mat']);
else
    infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
end
info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
% momentum file
momFile  = [info.rootMat,info.rootName,'MomentumTerms.nc'];
%
x   = ncread(momFile,'x');
y   = ncread(momFile,'y');
t   = ncread(momFile,'t');
%
if jj==1;
    PgradX = nan(length(x),Ndirs);
    RadStr = nan(length(x),Ndirs);
    Fbreak = nan(length(x),Ndirs);
    Xadvec = nan(length(x),Ndirs);
    Fdrag  = nan(length(x),Ndirs);
    %
    PgradX_rms = nan(length(x),Ndirs);
    RadStr_rms = nan(length(x),Ndirs);
    Fbreak_rms = nan(length(x),Ndirs);
    Xadvec_rms = nan(length(x),Ndirs);
    Fdrag_rms  = nan(length(x),Ndirs);
    Rx_rms  = nan(length(x),Ndirs);    
end
%
% load breaking force terms
BrkDissX = ncread(momFile,'BrkDissX');
DxSxx    = ncread(momFile,'DxSxx');
DySxy    = ncread(momFile,'DySxy');
dPdx     = ncread(momFile,'PgrdX');
DxUUH    = ncread(momFile,'DxUUH');
DyUVH    = ncread(momFile,'DyUVH');
FrcX     = ncread(momFile,'FRCX');
tmp      = dPdx+DxSxx+DySxy+DxUUH+DyUVH-BrkDissX;
%
BrkDissX_rms = std(BrkDissX,[],[1 3],'omitnan');
DxSxx_rms    = std(DxSxx   ,[],[1 3],'omitnan');
DySxy_rms    = std(DySxy   ,[],[1 3],'omitnan');
dPdx_rms     = std(dPdx    ,[],[1 3],'omitnan');
DxUUH_rms    = std(DxUUH   ,[],[1 3],'omitnan');
DyUVH_rms    = std(DyUVH   ,[],[1 3],'omitnan');
FrcX_rms     = std(FrcX    ,[],[1 3],'omitnan');
%
BrkDissX = mean(BrkDissX,[1 3],'omitnan');
DxSxx    = mean(DxSxx   ,[1 3],'omitnan');
DySxy    = mean(DySxy   ,[1 3],'omitnan');
dPdx     = mean(dPdx    ,[1 3],'omitnan');
DxUUH    = mean(DxUUH   ,[1 3],'omitnan');
DyUVH    = mean(DyUVH   ,[1 3],'omitnan');
FrcX     = mean(FrcX    ,[1 3],'omitnan');
tmp      = std( tmp     ,[],[1 3],'omitnan');
%
PgradX(:,jj) = dPdx;
RadStr(:,jj) = DxSxx+DySxy;
Fbreak(:,jj) = BrkDissX;
Xadvec(:,jj) = DxUUH+DyUVH;
Fdrag (:,jj) = FrcX;
PgradX_rms(:,jj) = dPdx_rms;
RadStr_rms(:,jj) = DxSxx_rms+DySxy_rms;
Fbreak_rms(:,jj) = BrkDissX_rms;
Xadvec_rms(:,jj) = DxUUH_rms+DyUVH_rms;
Fdrag_rms (:,jj) = FrcX_rms;
Rx_rms    (:,jj) = tmp;
T_INTV_mean(jj)   = info.T_INTV_mean;
%
end
%
% $$$ idx = find(x>50 & x<400);
Rx  = PgradX+RadStr+Xadvec+Fdrag-Fbreak;
% $$$ rmsErr   = sqrt( mean(Rx(idx,:).^2,1,'omitnan') )
% $$$ rmsPgrdX = sqrt( mean(PgradX(idx,:).^2,1,'omitnan') )
% create colormap
clrs = cmocean('thermal',Ndirs);
%
%
xm = 2.5;
ym = 2.5;
pw = 8;
ph = 2.5;
ag = 0.5;
ppos1 = [xm ym            pw ph];
ppos2 = [xm ym+ph+ag      pw ph];
ppos3 = [xm ym+2*(ph+ag)  pw ph];
ppos4 = [xm ym+3*(ph+ag)  pw ph];
ppos5 = [xm ym+4*(ph+ag)  pw ph];
cbpos = [xm+pw+ag ym ag ph*2/3];
ps    = [2*xm+pw+3*ag, 3*ym+5*ph+3*ag];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
%
colororder(clrs)
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x, Rx*1e3 ,'-','linewidth',2);
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$R_x$ (m/s)$^{2}\times 10^{-3}$','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','xlim',[50 400])
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x, -Fbreak*1e3 ,'-','linewidth',2);
xline(50,'--r')
ylabel('$-F_{\mathrm{br},x}$','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x, RadStr*1e3 ,'-','linewidth',2);
xline(50,'--r')
ylabel('$\partial_x S_{x,x}+\partial_y S_{x,y}$','interpreter','latex')
set(a3,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a4 = axes('units','centimeters','position',ppos4);
p4 = plot(x, Xadvec*1e3 ,'-','linewidth',2);
xline(50,'--r')
ylabel('$\partial_x UUH + \partial_y UVH$','interpreter','latex')
set(a4,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a5 = axes('units','centimeters','position',ppos5);
p5 = plot(x, PgradX*1e3 ,'-','linewidth',2);
xline(50,'--r')
ylabel('$gH\partial_x \eta$','interpreter','latex')
set(a5,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
cb = axes('units','centimeters','position',cbpos);
imagesc(0,[1:Ndirs],reshape(clrs,Ndirs,1,3))
xlabel('[s]','interpreter','latex')
set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
       'ticklabelinterpreter','latex','ytick',[1:Ndirs],'yticklabel',cellstr(num2str(T_INTV_mean')),'fontsize',12)
%
% $$$ a1.YAxis.Exponent=1e-2;
% $$$ a1.YAxis.Limits  = [-5 5]*1e-2;
% $$$ a1.YAxis.TickValues  = [-4:4]*1e-2;
% $$$ 
% $$$ a2.YAxis.Exponent=1e-2;
% $$$ a2.YAxis.Limits  = [-5 5]*1e-2;
% $$$ a2.YAxis.TickValues  = [-4:4]*1e-2;
% $$$ 
% $$$ a3.YAxis.Exponent=1e-2;
% $$$ a3.YAxis.Limits  = [-5 5]*1e-2;
% $$$ a3.YAxis.TickValues  = [-4:4]*1e-2;
% $$$ 
% $$$ a4.YAxis.Exponent=1e-2;
% $$$ a4.YAxis.Limits  = [-5 5]*1e-2;
% $$$ a4.YAxis.TickValues  = [-4:4]*1e-2;
% $$$ 
% $$$ a5.YAxis.Exponent=1e-2;
% $$$ a5.YAxis.Limits  = [-5 5]*1e-2;
% $$$ a5.YAxis.TickValues  = [-4:4]*1e-2;
% $$$ 
figname = [figDIR,'dominant_crossshore_momentum_budget_vs_TINTVmean.pdf'];
exportgraphics(fig,figname)
close(fig)
%
%
%
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
%
colororder(clrs)
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x, Rx_rms*1e2 ,'-','linewidth',2);
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('rms$(R_x)$ (m/s)$^{2}\times 10^{-2}$','interpreter','latex','fontsize',10)
set(a1,'tickdir','out','ticklabelinterpreter','latex','xlim',[50 400])
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x, Fbreak_rms*1e2 ,'-','linewidth',2);
xline(50,'--r')
ylabel('rms$(F_{\mathrm{br},x})$','interpreter','latex','fontsize',10)
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x, RadStr_rms*1e2 ,'-','linewidth',2);
xline(50,'--r')
ylabel('rms$(\partial_x S_{x,x}+\partial_y S_{x,y})$','interpreter','latex','fontsize',10)
set(a3,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a4 = axes('units','centimeters','position',ppos4);
p4 = plot(x, Xadvec_rms*1e2 ,'-','linewidth',2);
xline(50,'--r')
ylabel('rms$(\partial_x UUH + \partial_y UVH)$','interpreter','latex','fontsize',10)
set(a4,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a5 = axes('units','centimeters','position',ppos5);
p5 = plot(x, PgradX_rms*1e3 ,'-','linewidth',2);
xline(50,'--r')
ylabel('rms$(gH\partial_x \eta)$','interpreter','latex','fontsize',10)
set(a5,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
cb = axes('units','centimeters','position',cbpos);
imagesc(0,[1:Ndirs],reshape(clrs,Ndirs,1,3))
xlabel('[s]','interpreter','latex')
set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
       'ticklabelinterpreter','latex','ytick',[1:Ndirs],'yticklabel',cellstr(num2str(T_INTV_mean')),'fontsize',12)
figname = [figDIR,'dominant_crossshore_momentum_rms_vs_TINTVmean.pdf'];
exportgraphics(fig,figname)
close(fig)
