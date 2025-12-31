function fout = plot_FUNWAVE_run_Wave_statistics(info);
%
% USAGE: figure_list = plot_FUNWAVE_run_statistics(info);
%
% takes run "info" structure and creates plots of run
% statistics, outputing a list of archived figures.

figDIR = [info.rootMOD,filesep,'figures/'];
if ~exist(figDIR,'dir')
    eval(['!mkdir -p ',figDIR])
end
fout = {};
% 1) plot wave breaking statistics:
% 1.1) from info.waveForceFile:
%      Lc_pdf: (L,x) surface and line-plots
%      Nc: number of crests vs x (combine with L(x))
x   = ncread(info.waveForceFile,'x');
y   = ncread(info.waveForceFile,'y');
%
Xb     = ncread(info.waveForceFile,'Xb');
Lb     = ncread(info.waveForceFile,'Lb');
Lc_pdf = ncread(info.waveForceFile,'Lc_pdf');
Lc     = ncread(info.waveForceFile,'Lc');
Lc_std = ncread(info.waveForceFile,'Lc_plus_std');
Nc     = ncread(info.waveForceFile,'Lc');
%
%%  iii) line plots
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.25;
ppos1 = [xm ym           pw ph];
ppos2 = [xm ym+ph+ag     pw ph];
ppos3 = [xm 1.5*ym+2*ph+2*ag pw ph];
ps    = [2*xm+pw 3*ym+2*ag+3*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
plot(Xb,Nc,'-k','linewidth',2)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$N_\mathrm{c}$ [n/a]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
plot(Xb, Lc,'-k',Xb,Lc_std,'--r', 'linewidth',2)
ylabel('$L_\mathrm{c}$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
a3 = axes('units','centimeters','position',ppos3);
N  = floor(length(Xb)/10);
clrs = cmocean('thermal',N);
for jj=1:N
    tmp = mean(Lc_pdf(:,10*(jj-1) + [1:10]),2,'omitnan');
    hold on,
    loglog(Lb,tmp,'-','color',clrs(jj,:),'linewidth',2)
end
a3.XAxis.Scale='log';
a3.YAxis.Scale='log';
xlabel('$L$ [m]','interpreter','latex')
ylabel('$\mathrm{pdf}(L)$ [n/a]','interpreter','latex')
set(a3,'tickdir','out','ticklabelinterpreter','latex')
caxis([Xb(1) Xb(10*N)]),colormap(clrs),cb=colorbar;ylabel(cb,'$x$ [m]','interpreter','latex')
figname = [figDIR,info.runName,'_crest_length.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
% 1.2) from info.waveStats:
%      Hs: line vs x
%      Snn: alongshore average at sevaral x bins
Hs  = ncread(info.waveStatsFile,'Hs');
frq = ncread(info.waveStatsFile,'frq');
Snn = ncread(info.waveStatsFile,'Snn');
%
%%  iii) line plots
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.25;
ppos1 = [xm ym           pw ph];
ppos3 = [xm 2*ym+1*ph+1*ag pw ph];
ps    = [2*xm+pw 3*ym+1*ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
plot(x,mean(Hs,1,'omitnan'),'-k','linewidth',2)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$H_\mathrm{s}$ [m]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos3);
N  = floor(length(x)/20);
clrs = cmocean('thermal',N);
for jj=1:N
    tmp = mean(Snn(:,:,20*(jj-1) + [1:20]),[2 3],'omitnan');
    hold on,
    semilogx(a2,frq,tmp,'-','color',clrs(jj,:),'linewidth',2)
end
a2.XAxis.Scale='log';
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$S_{\eta \eta}$ [m$^2$]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex')
caxis(a2,[x(1) x(20*N)]), colormap(clrs),cb=colorbar;ylabel(cb,'$x$ [m]','interpreter','latex')
figname = [figDIR,info.runName,'_wave_height_and_spectra.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%
