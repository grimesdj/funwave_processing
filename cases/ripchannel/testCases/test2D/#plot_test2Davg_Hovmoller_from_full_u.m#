% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHY = 'test2D'
%
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
matDIR   = [runDIR,filesep,'mat_data'];
figDIR   = [runDIR,filesep,'figures/'];
%
% the list of run directories are saved in:
load([matDIR,filesep,'runs_to_process.mat'])
%
% filter for time averaging:
nf = 301;
flt=hanning(nf); flt=flt./sum(flt);
%
% loop over run_dirs
Ndirs  = length(run_dirs);
for jj = 1:Ndirs
% 1) get current run subdirectory to process:
runID    = run_dirs{jj};
fprintf('\n processing: %s %s \n', runBATHY,runID)    
% 2) get the archived info structure:
infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
%
% location of .nc files:
momFile  = [info.rootMOD,'mat_data',filesep,info.rootName,'MomentumTerms.nc'];
velFiles = dir([info.rootMOD,'mat_data',filesep,info.rootName,'u_*.nc']);
%
% estimate the location of the surfzone edge
x   = ncread([velFiles(1).folder,filesep,velFiles(1).name],'x');
y   = ncread([velFiles(1).folder,filesep,velFiles(1).name],'y');
% load breaking to estimate surfzone width
BrkDissX = ncread(momFile,'BrkDissX');
tmp0     = BrkDissX; rms0 = rms(tmp0,[1 3],'omitnan');
cum_rms  = cumsum(rms0);
idx_Lsz  = find(cum_rms>=0.98*max(cum_rms),1,'first');
Lsz(jj)  = x(idx_Lsz)
info.Lsz = x(idx_Lsz);
%
Nf       = length(velFiles);
t        = [];
HOV      = [];
for kk = 1:Nf
    velFile = [velFiles(kk).folder,filesep,velFiles(kk).name];
    y   = ncread(velFile,'y');
    t   = cat(1,t,ncread(velFile,'t'));
    tmp = ncread(velFile,'u',[1 idx_Lsz 1], [inf 100 inf]);
    HOV = cat(2,HOV,squeeze(mean(tmp,2)));
end
clear tmp
%
HOV = conv2(HOV,flt','same');
% estimate the alongshore wavenumber spectra of HOV
[Suu,ky] = alongshore_spectra_estimate(info,HOV);
Suu = mean(Suu,2);
Lrip= sum(Suu)./sum(Suu.*ky)
Sex(:,jj)= Suu;
%
% make a hovmoller plot
xm = 2.5;
ym = 2.5;
pw = 8;
ph = 3;
ag = 0.5;
ppos1 = [xm ym        pw ph];
ppos2 = [xm ym+ph+ym  pw ph];
cbpos = [xm+pw+ag ym ag ph*2/3];
ps    = [2*xm+pw+3*ag, 3.5*ym+2*ph+ag];
fig   = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
clims = 1.5*std(HOV(:))*[-1 1];
cm    = cmocean('balance');
clrs  = clims(1):diff(clims)/255:clims(2);
%
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,(t-t(1)),HOV')
colormap(cm), caxis(clims)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$t$ [s]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
%
a2 = axes('units','centimeters','position',ppos2);
semilogx(ky,Suu,'-k','linewidth',2)
xline(1/Lrip,'--r','linewidth',2)
xlabel('$k_y$ [m$^-1$]','interpreter','latex')
ylabel('$S_{u_\psi u_\psi}$ (m/s)$^2$','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex')
%
cb = axes('units','centimeters','position',cbpos);
imagesc(0,clrs,reshape(cm,256,1,3))
xlabel('[m/s]','interpreter','latex')
set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
       'ticklabelinterpreter','latex')
%
figname = [figDIR,info.runName,'_U_Hovmoller_and_Spectra.pdf'];
exportgraphics(fig,figname)
close(fig)
%
end
%
