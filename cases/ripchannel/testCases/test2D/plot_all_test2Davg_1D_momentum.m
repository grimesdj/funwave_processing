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
run_dirs = cat(1,{'plnr2D_h10t10s10d00'},run_dirs);
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
end
%
% load breaking force terms
BrkDissX = ncread(momFile,'BrkDissX');
DxSxx    = ncread(momFile,'DxSxx');
dPdx     = ncread(momFile,'PgrdX');
DxUUH    = ncread(momFile,'DxUUH');
FrcX     = ncread(momFile,'FRCX');
%
%
BrkDissX = mean(BrkDissX,[1 3],'omitnan');
DxSxx    = mean(DxSxx   ,[1 3],'omitnan');
dPdx     = mean(dPdx    ,[1 3],'omitnan');
DxUUH    = mean(DxUUH   ,[1 3],'omitnan');
FrcX     = mean(FrcX    ,[1 3],'omitnan');
%
PgradX(:,jj) = dPdx;
RadStr(:,jj) = DxSxx;
Fbreak(:,jj) = BrkDissX;
Xadvec(:,jj) = DxUUH;
Fdrag (:,jj) = FrcX;
T_INTV_mean(jj)   = info.T_INTV_mean;
%
end
%
idx = find(x>50 & x<400);
Rx  = PgradX+RadStr+Xadvec+Fdrag-Fbreak;
rmsErr   = sqrt( mean(Rx(idx,:).^2,1,'omitnan') )
rmsPgrdX = sqrt( mean(PgradX(idx,:).^2,1,'omitnan') )
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
cbpos = [xm+pw+ag ym ag ph*2/3];
ps    = [2*xm+pw+3*ag, 3*ym+4*ph+3*ag];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
%
colororder(clrs)
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x, Rx ,'-','linewidth',2);
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$R_x$ (m/s)$^{2}$','interpreter','latex')
a1.YAxis.Exponent = 1;
set(a1,'tickdir','out','ticklabelinterpreter','latex','xlim',[50 400])
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(x, -Fbreak ,'-','linewidth',2);
xline(50,'--r')
ylabel('$-F_{\mathrm{br},x}$','interpreter','latex')
a2.YAxis.Exponent = 1;
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x, RadStr ,'-','linewidth',2);
xline(50,'--r')
ylabel('$\partial_x S_{x,x}$','interpreter','latex')
a3.YAxis.Exponent = 1;
set(a3,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
a4 = axes('units','centimeters','position',ppos4);
p4 = plot(x, PgradX ,'-','linewidth',2);
xline(50,'--r')
ylabel('$gH\partial_x \eta$','interpreter','latex')
a4.YAxis.Exponent = 1;
set(a4,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',[50 400])
%
cb = axes('units','centimeters','position',cbpos);
imagesc(0,[1:Ndirs],reshape(clrs,Ndirs,1,3))
xlabel('[s]','interpreter','latex')
set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
       'ticklabelinterpreter','latex','ytick',[1:Ndirs],'yticklabel',cellstr(num2str(T_INTV_mean')),'fontsize',12)
%
figname = [figDIR,'dominant_crossshore_momentum_budget_vs_TINTVmean.pdf'];
exportgraphics(fig,figname)
close(fig)
