function fout = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info);
%
% USAGE: figure_list = plot_FUNWAVE_run_WaveAvgVelocity_statistics(info);
%
% takes run "info" structure and creates plots of run
% statistics, outputing a list of archived figures.
%
% also add stats to the velocity decomposition file:
% Uex: (mean) and (eddy) (m/s)
% Umke: mean kinetic energy (m/s)
% Ueke: eddy kinetic energy (m/s)
% ENSmean: mean enstrophy (1/s)^2
% ENSeddy: eddy enstrophy (1/s)^2

figDIR = [info.rootMOD,filesep,'figures/'];
if ~exist(figDIR,'dir')
    eval(['!mkdir -p ',figDIR])
end
fout = {};
%
fin = [info.rootMat,info.rootName,'dep.nc'];
h = ncread(fin,'dep');
%
%% 2) plot velocity and vorticity statistics
x        = ncread(info.rotVelFile,'x');
y        = ncread(info.rotVelFile,'y');
t        = ncread(info.rotVelFile,'t');
Urot     = ncread(info.rotVelFile,'Urot');
Vrot     = ncread(info.rotVelFile,'Vrot');
VORT     = ncread(info.rotVelFile,'VORT');
ETA      = ncread(info.rotVelFile,'eta');
H        = max(h+ETA,0.1);
%
%% estimate eddy energy/exchange statistics
tmp      = sqrt( Urot.^2 + Vrot.^2 );
Ueke      = sum( tmp.*H , [1 3],'omitnan')./sum( H, [1 3],'omitnan');
tmp      = Urot; tmp(Urot<0)=nan;
tmp1     = H;    tmp1(H<0.1)=nan;
Uex_eddy      = sum( tmp.*tmp1 , [1 3],'omitnan')./sum( tmp1, [1 3], 'omitnan');clear tmp
%
%% estimate the energy/exchange for the mean fields:
momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
Umean   = ncread(momFile,'umean');
Vmean   = ncread(momFile,'vmean');
ETAmean = ncread(momFile,'etamean');
%
% remove stokes velocity from mean:
H = h+ETAmean;
T = mean(Umean.*H,1);
Ustokes = T./mean(H,1);
Umean   = Umean-Ustokes;
Umean   = mean(Umean,3,'omitnan');
Vmean   = mean(Vmean,3,'omitnan');
ETAmean = mean(ETAmean,3,'omitnan');
Hmean   = h+ETAmean;
%
%% estimate mean kinetic energy
Umke     = sum(sqrt( Umean.^2 + Vmean.^2).*Hmean, 1,'omitnan')./sum(Hmean, 1,'omitnan');
%
tmp  = Umean; tmp(Umean<0)=nan;
tmp1 = Hmean; tmp1(Hmean<0.1)=nan;
Uex_mean = sum(tmp.*tmp1, 1, 'omitnan')./sum(tmp1, 1, 'omitnan');
%
%% estimate the enstrophy stats
[~,dUdy] = gradient(Umean./info.dy);
[dVdx,~] = gradient(Vmean./info.dx);
VORTavg  = dVdx-dUdy;
%
ENS_mean = sqrt( mean(VORTavg.^2, [1 3],'omitnan') );
ENS_eddy = sqrt( mean(VORT   .^2, [1 3],'omitnan') );
%
% figure properties
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
plot(x, sqrt(ENS_eddy),'-k',x, sqrt(ENS_mean),'--k', 'linewidth',2)
legend({'$\bar{\omega}$','$\langle\omega\rangle$'},'interpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\mathrm{rms}\langle \omega \rangle$ (s$^{-1}$)','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
a2 = axes('units','centimeters','position',ppos2);
plot(x, Uex_eddy,'-k',x,Uex_mean,'--k',x,Ueke,'-r',x,Umke,'-b', 'linewidth',2)
legend({'$U_\mathrm{ex}$','$\langle U\rangle_\mathrm{ex}$','$U_\mathrm{eke}$','$U_\mathrm{mke}$'},'interpreter','latex')
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
if ~isfield(info,'x_shoreline');
    alims = [x(info.subDomain(3:4))' y(info.subDomain(1:2))'];
else
    alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
end
%
%
scale = 1e-1;
if std(VORT(:))<scale
    rng = ceil(log10(3*std(VORT(:))/2));
    scale = 10^rng;
    clims = [-0.5 0.5]*scale;    
else
    clims = [-1 1]*1e-1;
end

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
    hold(ax0,'on'), contour(ax0,y,x,h',[0:1:6],'-k')
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\bar \\omega$ ',(mean(t(jj))-t(1))/60);
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
%% archive:
% check for existing vars: Uex_eddy, etc.
finfo = ncinfo(info.rotVelFile);
vars  = {finfo.Variables.Name}';
if ~ismember('Uex_eddy',vars)
dim_x  = {"x",length(x)};
nccreate  (info.rotVelFile,'Uex_eddy','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'Uex_eddy','Description','Exchange velocity from wave-averaged rotational velocity');
%
nccreate  (info.rotVelFile,'Uex_mean','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'Uex_mean','Description','Exchange velocity from time-mean velocity');
%
nccreate  (info.rotVelFile,'Ueke','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'Ueke','Description','Energy velocity scale from wave-averaged rotational velocity');
%
nccreate  (info.rotVelFile,'Umke','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'Umke','Description','Energy velocity scale from time-mean velocity');
%
nccreate  (info.rotVelFile,'ENS_eddy','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'ENS_eddy','Description','Enstrophy from wave-averaged rotational velocity');
%
nccreate  (info.rotVelFile,'ENS_mean','Dimensions',dim_x,'Format','netcdf4')
ncwriteatt(info.rotVelFile,'ENS_mean','Description','Enstrophy from time-mean velocity');
%
end
ncwrite   (info.rotVelFile,'Uex_eddy',Uex_eddy);
ncwrite   (info.rotVelFile,'Uex_mean',Uex_mean);
ncwrite   (info.rotVelFile,'Ueke',Ueke);
ncwrite   (info.rotVelFile,'Umke',Umke);
ncwrite   (info.rotVelFile,'ENS_eddy',ENS_eddy);
ncwrite   (info.rotVelFile,'ENS_mean',ENS_mean);
