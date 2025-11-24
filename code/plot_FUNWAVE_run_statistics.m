function fout = plot_FUNWAVE_run_statistics(info);
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
%      Ebr & Ebr_lp: line vs (x) and movie vs time
%      Ibr & Ibr_lp: same
%      Lc_pdf: (L,x) surface and line-plots
%      Nc: number of crests vs x (combine with L(x))
% 1.1a) start with Ebr and Ebr_lp
x   = ncread(info.waveForceFile,'x');
y   = ncread(info.waveForceFile,'y');
t   = ncread(info.waveForceFile,'t');
Ebr    = ncread(info.waveForceFile,'Ebr');
Ebr_lp = ncread(info.waveForceFile,'Ebr_lp');
%
%%    i) line plots
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.25;
ppos1 = [xm ym       pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [2*xm+pw 2*ym+ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
plot(x, mean(Ebr_lp,[1 3],'omitnan'), 'linewidth',2)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\langle \bar{u} \cdot F_\mathrm{br}\rangle$ ','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
plot(x, mean(Ebr,[1 3],'omitnan'), 'linewidth',2)
ylabel('$\langle u'' \cdot F_\mathrm{br}\rangle$ (m$^2$/s$^3$)','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
a1.YAxis.Exponent=0;
a2.YAxis.Exponent=0;
figname = [figDIR,info.runName,'_dissipation_rate.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%%    ii) make a movie
% define limits 
alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
clims = [-5 5]*1e-1;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [figDIR,info.runName,'_dissipation_rate_fast'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(Ebr(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle u'' \\cdot F_\\mathrm{br}\\rangle$ ',(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\langle \tilde{\epsilon}\,\rangle$ (m$^2$/s$^3$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
%
%
%
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [figDIR,info.runName,'_dissipation_rate_slow'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(Ebr_lp(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle \\bar{u} \\cdot F_\\mathrm{br}\\rangle$ ',(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\langle \bar{\epsilon}\,\rangle$ (m$^2$/s$^3$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
%
clear Ebr_lp Ebr
%
%% next, plot the enstrophy production rate
VProd    = ncread(info.waveForceFile,'VProd');
VProd_lp = ncread(info.waveForceFile,'VProd_lp');
%
ppos1 = [xm ym       pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [2*xm+pw 2*ym+ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
plot(x, mean(VProd_lp,[1 3],'omitnan'), 'linewidth',2)
xlabel('$x$ [m]','interpreter','latex','fontsize',10)
ylabel('$\langle \bar{\omega} \mathrm{curl}(F_\mathrm{br})\rangle$ ','interpreter','latex','fontsize',10)
a1.YAxis.Exponent=0;
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
plot(x, mean(VProd,[1 3],'omitnan'), 'linewidth',2)
ylabel('$\langle \omega'' \mathrm{curl}(F_\mathrm{br})\rangle$ (s$^{-3}$)','interpreter','latex')
a2.YAxis.Exponent=0;
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
figname = [figDIR,info.runName,'_enstrophy_production_rate.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%%    ii) make a movie
% define limits 
alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
clims = [-1 1]*1e-2;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [figDIR,info.runName,'_enstrophy_production_rate_fast'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(VProd(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle \\omega'' \\mathrm{curl}(F_\\mathrm{br})\\rangle$ ',(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\langle \tilde{\Omega}\,\rangle$ (s$^{-3}$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
%
%
%
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [figDIR,info.runName,'_enstrophy_production_rate_slow'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(VProd_lp(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle \\bar{\\omega} \\mathrm{curl}(F_\\mathrm{br})\\rangle$ ',(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\langle \bar{\Omega}\,\rangle$ (s$^{-3}$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
%
clear VProd_lp VProd
%
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
%
ky      = ncread(info.waveForceFile,'ky');
Icoh    = ncread(info.waveForceFile,'Icoh');
Icoh_lp = ncread(info.waveForceFile,'Icoh_lp');
%
%%    i) line plots
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.25;
ppos1 = [xm ym       pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [2*xm+pw 2*ym+ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
N  = floor(length(x)/20);
clrs = cmocean('thermal',N);
for jj=1:N
    tmp = squeeze(mean(Icoh_lp(:,20*(jj-1) + [1:20],:),[2 3],'omitnan'));
    hold on,
    semilogx(a1,ky, tmp,'color',clrs(jj,:), 'linewidth',2)
end
a1.XAxis.Scale='log';
xlabel('$k_y$ [1/m]','interpreter','latex')
ylabel('$\mathrm{coh}(\omega'',F_\mathrm{br})$ [n/a]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
Icoh_lp=mean(Icoh_lp,3,'omitnan');
for jj=1:N
    tmp = squeeze(mean(Icoh_lp(:,20*(jj-1) + [[1:20]]),2,'omitnan'));
    hold on,
    semilogx(a2,ky, tmp,'color',clrs(jj,:), 'linewidth',2)
end
a2.XAxis.Scale='log';
xlabel('$k_y$ [1/m]','interpreter','latex')
ylabel('$\mathrm{coh}(\bar{\omega},F_\mathrm{br})$ [n/a]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
caxis(a2,[x(1) x(20*N)]), colormap(clrs),cb=colorbar;ylabel(cb,'$x$ [m]','interpreter','latex')
figname = [figDIR,info.runName,'_enstrophy_production_coherence.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
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
%% 2) plot velocity and vorticity statistics
Urot     = ncread(info.rotVelFile,'Urot');
Vrot     = ncread(info.rotVelFile,'Vrot');
VORT     = ncread(info.rotVelFile,'VORT');
%
% estimate energy/exchange statistics
tmp      = sqrt( Urot.^2 + Vrot.^2 );
EKE      = mean( tmp , [1 3],'omitnan');
tmp      = Urot; tmp(Urot<0)=nan;
Uex      = mean( tmp , [1 3],'omitnan');clear tmp
%
ppos1 = [xm ym       pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [2*xm+pw 2*ym+ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
a1 = axes('units','centimeters','position',ppos1);
plot(x, rms(VORT,[1 3],'omitnan'), 'linewidth',2)
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\mathrm{rms}\langle \omega \rangle$ (s$^{-1}$)','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
plot(x, Uex,'-k',x,EKE,'-r', 'linewidth',2)
legend({'$U_\mathrm{ex}$','$U_\mathrm{eke}$'},'interpreter','latex')
ylabel('(m/s)','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
figname = [figDIR,info.runName,'_velocity_and_vorticity_stats.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%%    ii) make a movie
% define limits 
alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
clims = [-1 1]*1e-1;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [figDIR,info.runName,'_vorticity'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(VORT(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle \\omega \\rangle$ ',(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'(s$^{-1}$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
%
%
clear VORT Urot Vrot EKE Urot
%
