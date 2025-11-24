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
    rmsFbr_visc = nan(length(x),Ndirs);
    rmsFbr_rad  = nan(length(x),Ndirs);
% $$$ elseif jj==Ndirs
% $$$     % breaking force files
% $$$     FbrXFiles  = dir([info.rootMat,info.rootName,'BrkSrcX_*.nc']);
% $$$     FbrYFiles  = dir([info.rootMat,info.rootName,'BrkSrcY_*.nc']);
% $$$     rmsFbr_visc0 = zeros(1,length(x));
% $$$     for kk=1:length(FbrXFiles)
% $$$         % load breaking force terms
% $$$         BrkDissX = ncread([FbrXFiles(kk).folder,filesep,FbrXFiles(kk).name],'BrkSrcX');
% $$$         BrkDissY = ncread([FbrYFiles(kk).folder,filesep,FbrYFiles(kk).name],'BrkSrcY');
% $$$         %
% $$$         [dyFx,~   ] = gradientDG(BrkDissX/info.dy);
% $$$         [~   ,dxFy] = gradientDG(BrkDissY/info.dx);
% $$$         cFbr        = dxFy - dyFx;
% $$$         rmsFbr_visc0 = rmsFbr_visc0 + rms(cFbr,[1 3],'omitnan').^2;
% $$$     end
% $$$     rmsFbr_visc0 = sqrt(rmsFbr_visc0/length(FbrXFiles));
end
%
% load breaking force terms
BrkDissX = ncread(momFile,'BrkDissX');
BrkDissY = ncread(momFile,'BrkDissY');
DxSxx    = ncread(momFile,'DxSxx');
DySyy    = ncread(momFile,'DySyy');
%
[dyFx,~   ] = gradientDG(BrkDissX/info.dy);
[~   ,dxFy] = gradientDG(BrkDissY/info.dx);
cFbr        = dxFy - dyFx;
[dyFx,~   ] = gradientDG(DxSxx/info.dy);
[~   ,dxFy] = gradientDG(DySyy/info.dx);
cS          = dxFy - dyFx;
%
tmp0 = cFbr; tmp0(~info.mask)=nan;  rms0 = rms(tmp0,[1 3],'omitnan');
tmp1 = cS;   tmp1(~info.mask)=nan;  rms1 = rms(tmp1,[1 3],'omitnan');
%
rmsFbr_visc(:,jj) = rms0;
rmsFbr_rad (:,jj) = rms1;
T_INTV_mean(jj)   = info.T_INTV_mean;
%
end
%
% create colormap
clrs = cmocean('thermal',Ndirs);
%
%
xm = 2.5;
ym = 2.5;
pw = 8;
ph = 3;
ag = 0.5;
ppos1 = [xm ym       pw ph];
cbpos = [xm+pw+ag ym ag ph*2/3];
ps    = [2*xm+pw+3*ag, 2.5*ym+ph];
fig   = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
colororder(clrs)
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,rmsFbr_visc,'-',x,rmsFbr_rad,'--','linewidth',2);
hold on,
p2 = plot(x,rmsFbr_visc0,':r','linewidth',2)
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('rms() [m/s$^{2}$]','interpreter','latex')
legend([p1([1 Ndirs+1]); p2],{'curl$(\bar{F}_\mathrm{br})$','curl$(\nabla S)$','curl$({F}_\mathrm{br})$'},'interpreter','latex','fontsize',10)
% legend([p1([1 Ndirs+1])],{'curl$(\bar{F}_\mathrm{br})$','curl$(\nabla S)$'},'interpreter','latex','fontsize',10)
% a1.YAxis.YScale = 'log';
set(a1,'tickdir','out','ticklabelinterpreter','latex')
%
cb = axes('units','centimeters','position',cbpos);
imagesc(0,[1:Ndirs],reshape(clrs,Ndirs,1,3))
xlabel('[s]','interpreter','latex')
set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
       'ticklabelinterpreter','latex','ytick',[1:Ndirs],'yticklabel',cellstr(num2str(T_INTV_mean')),'fontsize',12)
%
figname = [figDIR,'rms_curl_Fbr_all.pdf'];
exportgraphics(fig,figname)
close(fig)
