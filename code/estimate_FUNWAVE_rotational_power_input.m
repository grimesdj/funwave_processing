function info = estimate_FUNWAVE_rotational_power_input(info)
%
% USAGE: info = estimate_FUNWAVE_rotational_power_input(info)
%
disp(['working on : ', info.runName])
rot_vel = matfile(info.rotVelFile);
rot_fbr = matfile(info.waveForceDecompFile);
force   = matfile(info.waveForceFile);
x       = force.x; dx = x(2)-x(1);
y       = force.y; dy = y(2)-y(1);
t       = force.t; t = t-t(1);
Nx      = length(x);
Ny      = length(y);
Nt      = length(t);
dt      = mean(diff(t));
% for mapping to the local shoreline coordinates:
xsl        = info.x_shoreline';
xp         = x-mean(xsl);
x_map      = x-xsl;
[xx,yy]    = meshgrid(x-mean(xsl),y);
%
%
% 0) create file for saving stats
power = matfile([info.rootMat,info.rootName,'rotational_power_input.mat'],'writable',true);
info.rotPowerFile = power.Properties.Source;
save(info.fileName,'-struct','info')
%
% 1) estimate Prot = Urot * Fbr
Prot = force.Fx.*rot_vel.Urot + force.Fy.*rot_vel.Vrot;
% $$$ Prot = power.Prot;
%
% 2) estimate Prot_avg = < Prot > = <Urot>*<Fbr> + <(Urot-<Urot>)*(Fbr-<Fbr>)>
Prot_avg = mean(Prot,3);
Urot_avg = mean(rot_vel.Urot,3);
Vrot_avg = mean(rot_vel.Vrot,3);
Fx_avg   = mean(force.Fx,3);
Fy_avg   = mean(force.Fy,3);
Prot_avg1= Urot_avg.*Fx_avg + Vrot_avg.*Fy_avg;
Prot_avg2= mean( (rot_vel.Urot-Urot_avg).*(force.Fx-Fx_avg) + (rot_vel.Vrot-Vrot_avg).*(force.Fy-Fy_avg), 3);
% 3) estimate Prot_var = Prot - Prot_avg = <Urot>*(Fbr-<Fbr>) + (Urot-<Urot>)*<Fbr> + ((Urot-<Urot>)(Fbr-<Fbr>) - <(Urot-<Urot>)(Fbr-<Fbr>)>)
Prot_var = Prot - Prot_avg;
Prot_var1= (Urot_avg      .* force.Fx  + Vrot_avg     .* force.Fy) - Prot_avg1;
Prot_var2= (rot_vel.Urot  .* Fx_avg    + rot_vel.Vrot .* Fy_avg  ) - Prot_avg1;
% log basic output
power.x    = x;
power.y    = y;
power.t    = t;
power.Prot      = Prot;
power.Prot_avg  = Prot_avg;
power.Prot_avg1 = Prot_avg1;
power.Prot_avg2 = Prot_avg2;
power.Prot_var  = Prot_var;
power.Prot_var1 = Prot_var1;
power.Prot_var2 = Prot_var2;
%
% 4a) make 3-panel map of Prot_avg and components,
xm = 2; ym = 2; pw = 10; ph = 3; ag=0.25; ps = [1.5*xm+pw, 2*ym+4*ag+3*ph];
ppos1 = [xm ym             pw ph];
ppos2 = [xm ym+ph+ag       pw ph];
ppos3 = [xm ym+2*ph+2*ag   pw ph];
cbpos = [xm ym+3*ph+3*ag 0.25*pw ag/2];
fig0  = figure('units','centimeters');
pos   = get(fig0,'position');
pos(3:4) = ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
clims= 0.1*[-1 1];
cm   = cmocean('balance');
clrs = clims(1):range(clims)/255:clims(2);
% mean field
f0ax1  = axes('units','centimeters','position',ppos3);
imagesc(f0ax1,y',xp',Prot_avg'), caxis(f0ax1,clims), colormap(f0ax1,cm),
hold(f0ax1,'on'), plot(f0ax1,y,xsl-mean(xsl),'--k')
ylabel(f0ax1,'$x$ [m]','interpreter','latex')
text(f0ax1,1250,145,'c) $\langle u_\mathrm{rot} \cdot F_\mathrm{br}\rangle$','interpreter','latex','fontsize',12)
title(info.runName)
set(f0ax1,'tickdir','out','ticklabelinterpreter','latex','ylim',-25+[0 200],'ydir','normal','xdir','reverse','xticklabel',[])
% product of means
f0ax2  = axes('units','centimeters','position',ppos2);
imagesc(f0ax2,y',xp',Prot_avg1'), caxis(f0ax2,clims), colormap(f0ax2,cm),
hold(f0ax2,'on'), plot(f0ax2,y,xsl-mean(xsl),'--k')
ylabel(f0ax2,'$x$ [m]','interpreter','latex')
text(f0ax2,1250,145,'b) $\langle u_\mathrm{rot}\rangle \cdot \langle F_\mathrm{br}\rangle$','interpreter','latex','fontsize',12)
set(f0ax2,'tickdir','out','ticklabelinterpreter','latex','ylim',-25+[0 200],'ydir','normal','xticklabel',[],'xdir','reverse')
% product of fluctuations
f0ax3  = axes('units','centimeters','position',ppos1);
imagesc(f0ax3,y',xp',Prot_avg2'), caxis(f0ax3,clims), colormap(f0ax3,cm),
hold(f0ax3,'on'), plot(f0ax3,y,xsl-mean(xsl),'--k')
ylabel(f0ax3,'$x$ [m]','interpreter','latex')
xlabel(f0ax3,'$y$ [m]','interpreter','latex')
text(f0ax3,1250,145,'a) $ \langle \tilde{u}_\mathrm{rot} \cdot \tilde{F}_\mathrm{br}\rangle$','interpreter','latex','fontsize',12)
set(f0ax3,'tickdir','out','ticklabelinterpreter','latex','ylim',-25+[0 200],'ydir','normal','xdir','reverse')
% make colorbar
f0cb = axes('units','centimeters','position',cbpos);
imagesc(f0cb,clrs,0,reshape(cm,1,256,3))
xlabel(f0cb,'[W/kg]','interpreter','latex','rotation',0,'fontsize',12)
set(f0cb, 'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','top','ytick',[],'ticklength',2*get(f0cb,'ticklength'))
% save figure
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_time_mean_vs_xy.pdf'])
%
% 4b) plot y-avg[], y-rms[], and y-spectra[] of Prot_avg
% first, map to shoreline coordinates
if std(xsl)>=3
    fprintf('\tmapping to shoreline coordinates\n')
    Prot_avg_yavg = griddata(x_map,yy,Prot_avg,xx,yy);
else
    Prot_avg_yavg = Prot_avg;
end
% estimate spectra
[Prot_avg_spec,ky] = alongshore_spectra_estimate(info,Prot_avg_yavg);
% estimate y-rms
Prot_avg_yrms = rms(Prot_avg_yavg,1,'omitnan');
% estimate y-avg
Prot_avg_yavg = nanmean(Prot_avg_yavg,1);
%
power.Prot_avg_yavg = Prot_avg_yavg;
power.Prot_avg_yrms = Prot_avg_yrms;
power.Prot_avg_spec = Prot_avg_spec;
%
% now plot all three, binning the spectra into five cross-shore bins
xm = 2; ym = 2; pw = 10; ph = 3.5; ag=0.25; ps = [1.5*xm+pw, 3*ym+2*ag+2.5*ph];
ppos1 = [xm        ym             pw       1.5*ph];
ppos2 = [xm        2*ym+1.5*ph        pw       1*ph];
cbpos = [xm+0.6*pw ym+1.25*ph      0.3*pw   ag];
fig1  = figure('units','centimeters');
pos   = get(fig1,'position');
pos(3:4) = ps;
set(fig1,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
% y-avg & y-rms fields
f1ax1  = axes('units','centimeters','position',ppos2);
plot(f1ax1,xp, Prot_avg_yavg,'-k',xp,Prot_avg_yrms,'-r','linewidth',2)
xlabel(f1ax1,' $(x-\bar{x}_\mathrm{sl})$~[m]','interpreter','latex')
ylabel(f1ax1,' ~[W/kg]','interpreter','latex')
f1ax1.YAxis.Exponent=-2;
f1l1  = legend('$\overline{\langle u_\mathrm{rot}\cdot F_\mathrm{br} \rangle}$','$\overline{{\langle u_\mathrm{rot}\cdot F_\mathrm{br} \rangle}^2}^{1/2}$');
set(f1l1,'interpreter','latex','fontsize',12)
title(info.runName)
set(f1ax1,'tickdir','out','ticklabelinterpreter','latex','xlim',-25+[0 300])
%
f1ax2  = axes('units','centimeters','position',ppos1);
hold(f1ax2,'on')
% now cross-shore average into 25m bins (dof~2*50, less the cross-shore correlation)
db    = 50;
xbins = [0:db:150]+db/2;
Nb    = length(xbins);
cm2   = cmocean('thermal',Nb);
Prot_avg_spec_binned = nan(length(ky),Nb);
for ii = 1:Nb
    inbin  =  find( xp>=(xbins(ii)-db/2) & xp<(xbins(ii)+db/2));
    tmp    = nanmean(Prot_avg_spec(:,inbin),2);
    loglog(f1ax2,ky,sqrt(tmp*(ky(2)-ky(1))),'color',cm2(ii,:),'linewidth',2),
    Prot_avg_spec_binned(:,ii)=tmp;
end
%
power.xbins = xbins;
power.ky    = ky;
power.Prot_avg_spec_binned = Prot_avg_spec_binned;
%
xlabel(f1ax2,' $k_y$~[m$^{-1}$]','interpreter','latex')
ylabel(f1ax2,' [W/kg]','interpreter','latex')
grid(f1ax2,'on')
set(f1ax2,'tickdir','out','ticklabelinterpreter','latex','box','on')
f1ax2.XAxis.Scale = 'log';
% $$$ f1ax2.YAxis.Scale = 'log';
%
f1cb = axes('units','centimeters','position',cbpos);
imagesc(f1cb,xbins,0,reshape(cm2,1,Nb,3))
xlabel(f1cb,'$(x-\bar{x}_\mathrm{sl})$','interpreter','latex','rotation',0,'fontsize',12)
set(f1cb, 'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','bottom','ytick',[],'ticklength',2*get(f1cb,'ticklength'),'xtick',xbins)
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_time_mean_alongshore_stats_vs_x.pdf'])
%
%
% 5a) make movie frame of wave-averaged (30-s or 60-s) Prot_var and estimate fractional contributions from each component,
% ***Only generate video object if we can make several 60s averages***
Navg   = round(min(60/dt,Nt));
Nplt   = round(Navg/2);
Nframes= floor(Nt/Nplt)-1;
% preallocate alongshore averaged stats
Prot_var_yavg = nan(1,Nx,Nframes);
Prot_var_yrms = nan(1,Nx,Nframes);
Prot_var_spec_binned = nan(length(ky),Nb,Nframes);
if Nframes>10
%
alims = [75 400 0 1500];
clims = [-1 1]*1e-2;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig2,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
% $$$ ps  = ps + [0 0.25];
% $$$ pos = pos+[0 0 0 0.25];
% $$$ set(fig2,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
% $$$ pos0 = get(ax0,'position');
% $$$ set(ax0,'position',pos0+[0 -0.75 0 0])
% $$$ pos00 = get(ax00,'position');
delete(ax00)
% $$$ pos01 = get(cx01,'position');
% $$$ set(cx01,'position',pos01+[0 -0.75 0 0])
%
%
vidName= sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_avg%ds'],round(Navg*dt));
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
mask = info.mask;
% cross-shore average to de-alias the under-resolved propagation
Nfx = round(10/dx); if ~mod(Nfx,2), Nfx = Nfx+1; end
FILT= hamming(Nfx); FILT = FILT/sum(FILT);
for jj=1:Nframes
    ii = Nplt*(jj-1) + [1:Navg];
    % estimate avg of total input, and fractional input from the mean, and component of varying components
    Prot_temp = Prot(:,:,ii);
    Prot_temp_avg = sum(nanmean(Prot_temp,3).*mask,[1 2],'omitnan');
    % fraction from mean
    frac_Prot_avg  = 1-sum(nanmean(Prot_temp-Prot_avg         ,3).*mask,[1 2],'omitnan')./Prot_temp_avg;
    % fraction from 1st varying term
    frac_Prot_var1 = 1-sum(nanmean(Prot_temp-Prot_var1(:,:,ii),3).*mask,[1 2],'omitnan')./Prot_temp_avg;
    % fraction from 2nd varying term    
    frac_Prot_var2 = 1-sum(nanmean(Prot_temp-Prot_var2(:,:,ii),3).*mask,[1 2],'omitnan')./Prot_temp_avg;
    %
    Prot_var_avg   = nanmean(Prot_var(:,:,ii),3);%nanmean(Eout(:,:,ii),3);
    Prot_var_smooth= conv2(Prot_var_avg,FILT','same');
    % plot avg
    imagesc(ax0,y,x,Prot_var_smooth'), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle u_\\mathrm{rot} \\cdot F_\\mathrm{br}\\rangle$ = %1.1f\\%%, $\\langle u_\\mathrm{rot}\\rangle \\cdot \\tilde{F}_\\mathrm{br}$ = %1.1f\\%%, $\\tilde{u}_\\mathrm{rot}\\cdot \\langle F_\\mathrm{br}\\rangle$ = %1.1f\\%% ',(mean(t(ii))-t(1))/60, 100*frac_Prot_avg, 100*frac_Prot_var1, 100*frac_Prot_var2);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\langle \tilde{\epsilon}\,\rangle_\mathrm{60 s}$ [W/kg]'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    %
    % get+write frame
    set(fig2,'position',pos);
    frame = getframe(fig2);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
    % now estimate alongshore stats on the wave-average
    % first, map to shoreline coordinates
    if std(xsl)>=3
        fprintf('\tmapping to shoreline coordinates\n')
        Prot_var_yavg_tmp = griddata(x_map,yy,Prot_var_smooth,xx,yy);
    else
        Prot_var_yavg_tmp = Prot_var_smooth;
    end
    % estimate spectra
    [Prot_var_spec_tmp,~] = alongshore_spectra_estimate(info,Prot_var_yavg_tmp);
    % estimate y-rms
    Prot_var_yrms_tmp = rms(Prot_var_yavg_tmp,1,'omitnan');
    % estimate y-avg
    Prot_var_yavg_tmp = nanmean(Prot_var_yavg_tmp,1);
    %
    Prot_var_yavg(1,:,jj) = Prot_var_yavg_tmp;
    Prot_var_yrms(1,:,jj) = Prot_var_yrms_tmp;
    %
    Prot_var_spec_tmp_binned = nan(length(ky),Nb);
    for kk = 1:Nb
        inbin  =  find( xp>=(xbins(kk)-db/2) & xp<(xbins(kk)+db/2));
        tmp    = nanmean(Prot_var_spec_tmp(:,inbin),2);
        Prot_var_spec_tmp_binned(:,kk)=tmp;
    end
    Prot_var_spec_binned(:,:,jj) = Prot_var_spec_tmp_binned;        
end


end
% 5b) plot y-avg[], y-rms[], and y-spectra[] of Prot_var
power.Prot_var_yavg = Prot_var_yavg;
power.Prot_var_yrms = Prot_var_yrms;
power.Prot_var_spec_binned = Prot_var_spec_binned;
%
%
% now plot all three, binning the spectra into five cross-shore bins
% now plot all three, binning the spectra into five cross-shore bins
xm = 2; ym = 2; pw = 10; ph = 3.5; ag=0.25; ps = [1.5*xm+pw, 3*ym+2*ag+2.5*ph];
ppos1 = [xm        ym             pw       1.5*ph];
ppos2 = [xm        2*ym+1.5*ph        pw       1*ph];
cbpos = [xm+0.6*pw ym+1.25*ph      0.3*pw   ag];
fig3  = figure('units','centimeters');
pos   = get(fig3,'position');
pos(3:4) = ps;
set(fig3,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
% y-avg & y-rms fields
f3ax1  = axes('units','centimeters','position',ppos2);
plot(f3ax1,xp, nanmean(Prot_var_yavg,3),'-k',xp,nanmean(Prot_var_yrms,3),'-r','linewidth',2)
xlabel(f3ax1,' $(x-\bar{x}_\mathrm{sl})$~[m]','interpreter','latex')
ylabel(f3ax1,' ~[W/kg]','interpreter','latex')
f3ax1.YAxis.Exponent=-2;
f3l1  = legend('$\overline{\langle \tilde{u}_\mathrm{rot}\cdot \tilde{F}_\mathrm{br} \rangle_\mathrm{60~s}}$','$\overline{{\langle \tilde{u}_\mathrm{rot}\cdot \tilde{F}_\mathrm{br} \rangle_\mathrm{60~s}}^2}^{1/2}$');
set(f3l1,'interpreter','latex','fontsize',12)
title(info.runName)
set(f3ax1,'tickdir','out','ticklabelinterpreter','latex','xlim',-25+[0 300])
%
f3ax2  = axes('units','centimeters','position',ppos1);
hold(f3ax2,'on')
% now cross-shore average into 25m bins (dof~2*50, less the cross-shore correlation)
for jj = 1:Nb
    tmp = nanmean(Prot_var_spec_binned(:,jj,:),3);
    loglog(f3ax2,ky,sqrt(tmp*(ky(2)-ky(1))),'color',cm2(jj,:),'linewidth',2),
end
%
xlabel(f3ax2,' $k_y$~[m$^{-1}$]','interpreter','latex')
ylabel(f3ax2,' [W/kg]','interpreter','latex')
grid(f3ax2,'on')
set(f3ax2,'tickdir','out','ticklabelinterpreter','latex','box','on')
f3ax2.XAxis.Scale = 'log';
% $$$ f3ax2.YAxis.Scale = 'log';
%
f3cb = axes('units','centimeters','position',cbpos);
imagesc(f3cb,xbins,0,reshape(cm2,1,Nb,3))
xlabel(f3cb,'$(x-\bar{x}_\mathrm{sl})$','interpreter','latex','rotation',0,'fontsize',12)
set(f3cb, 'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','xaxislocation','bottom','ytick',[],'ticklength',2*get(f3cb,'ticklength'),'xtick',xbins)
exportgraphics(fig3,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_time_varying_alongshore_stats_vs_x.pdf'])
%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now, what statistics to analyze?
[Ny,Nx,Nt] = size(Prot);
[nx,ny]    = meshgrid(1:Nx,1:Ny);
%
% map forcing to distance along-crest [0:0.01:1]?
l          = [0:0.01:1];
Eout       = nan(Ny,Nx,Nt);
for jj=1:Nt
    E  = Prot(:,:,jj);
    BW = bwboundaries(abs(E)>0.001);
    Nb = length(BW);
    for nb = 1:Nb
        points = BW{nb};
        yb     = points(:,1);
        xb     = points(:,2);
        % get grid that's inside current front
        in     = inpolygon(nx,ny,xb,yb);        
        % points inside polygon
        Ii     = find(in);
        xi     = nx(in);
        yi     = ny(in);
        ei     = E(in);
        % get the x-average magnitude and direction
        % first, get unique set of y-coords
        yu     = unique(yi);
        Nu     = length(yu);
        for nu = 1:Nu
            iu = find(yi==yu(nu));
            e_tot    = nanmean(ei(iu));
            ei(iu)   = e_tot;
            % ei(iu)   = real(e_tot)+imag(e_tot);
        end
        E(in) = ei;
% $$$         % get the length of the front
% $$$         lb(nb) = range(yu);
% $$$         mag(nb,:) = interp(0:1/(Nu-1):1,mag_uFbr,l);
% $$$         uFbr(nb,:) = interp(0:1/(Nu-1):1,mag_uFbr,l);
    end
    Eout(:,:,jj)=real(E);
end
%

%
% generate video object? only if we can make several 60s averages
Navg   = round(60/dt);
Nplt   = round(Navg/4);
Nframes= floor(Nplt*Nt/Navg);
if Nframes>10
%
alims = [75 400 0 1500];
clims = [-1 1]*1e-2;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig0,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
ps  = ps + [0 0.25];
pos = pos+[0 0 0 0.25];
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
pos0 = get(ax0,'position');
set(ax0,'position',pos0+[0 -0.75 0 0])
pos00 = get(ax00,'position');
set(ax00,'position',pos00+[0 -0.25 3 1.5])
pos01 = get(cx01,'position');
set(cx01,'position',pos01+[0 -0.75 0 0])
%
%
vidName= sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_avg%ds'],round(Navg*dt));
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
% E    = Eout(:,:,1);
Enz  = (Prot~=0);
% Etot = squeeze(nansum(nansum(Eout,1),2)./sum(sum(Enz,1),2));clear Enz
Etot = squeeze(nansum(nansum(Prot,1),2)./sum(sum(Enz,1),2));clear Enz
% smooth Etot using 20-second running mean
Etot_avg = conv(Etot,hamming(Navg)/sum(hamming(Navg)),'same');
for jj=1:Navg/Nplt:Nframes
    ii = Navg/Nplt*(jj-1) + [1:Navg];
    E  = nanmean(Prot(:,:,ii),3);%nanmean(Eout(:,:,ii),3);
    % plot avg
    imagesc(ax0,y,x,E'), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'$\epsilon=\left<u_\mathrm{rot}\cdot F_\mathrm{br}\right>$ [m$^2$/s$^{3}$]'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % plot the time evolution of
    yline(ax00,0,'--k')
    plot(ax00,t/60,Etot,'k.',t/60,Etot_avg,'-b',t(jj+Navg/2)/60,Etot_avg(jj+Navg/2),'ro','markersize',10);
    set(ax00,'xlim',(t(jj+Navg/2)+1*Navg*[-1 1])/60,'tickdir','out','ylim',[min(Etot_avg) max(Etot_avg)])
    xlabel(ax00,'$t$ [min]','interpreter','latex')
    ylabel(ax00,'$\overline{\epsilon}$ [m$^2$/s$^3$]','interpreter','latex')
    %
    % get+write frame
    set(fig0,'position',pos);
    frame = getframe(fig0);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
end
end
%
% now plot the time-averaged domain rotational power input
Eavg  = nanmean(Prot,3);
Eavg_smooth= conv2(Eavg,FILT','same');
%
alims = [75 400 0 1500];
clims = [-1 1]*1e-2;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
xm = 2; ym = 2; pw = 10; ph = 3; ag=0.25; ps = [1.5*xm+pw, 2*ym+ag+ph];
ppos = [xm ym             pw ph];
cbpos = [xm ym+ph+ag 0.25*pw ag/1.5];
fig0  = figure('units','centimeters');
pos   = get(fig0,'position');
pos(3:4) = ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
%
ax0  = axes('units','centimeters','position',ppos);
imagesc(ax0,y,x-mean(info.x_shoreline),Eavg_smooth'),
hold(ax0,'on'),plot(ax0,y,info.x_shoreline-mean(info.x_shoreline),'--k','linewidth',1.5)
caxis(ax0,clims)
colormap(ax0,cm)
ylabel(ax0,'$x$ [m]','interpreter','latex')
xlabel(ax0,'$y$ [m]','interpreter','latex')
set(ax0,'tickdir','out','ticklabelinterpreter','latex','ylim',-25+[0 200],'ydir','normal','xdir','reverse','xtick',[0 500 1000])
%
% make colorbar
cx01 = axes('units','centimeters','position',cbpos);
imagesc(cx01,clrs,0,reshape(cm,1,256,3))
% ylabel(cx01,{sprintf('$~~\\left<\\vec{u}_\\mathrm{rot}\\cdot \\vec{F}_\\mathrm{br}\\right>_\\mathrm{%d~s}$ [m$^2$/s$^{3}$]',round(range(t)))},'interpreter','latex','rotation',0,'fontsize',12,'verticalalignment','bottom','horizontalalignment','left')
ylabel(cx01,{'$~~\left<\vec{u}_\mathrm{rot}\cdot \vec{F}_\mathrm{br}\right>$ [m$^2$/s$^{3}$]'},'interpreter','latex','rotation',0,'fontsize',12,'verticalalignment','bottom','horizontalalignment','left')
set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','yaxislocation','right')
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_time_averaged.pdf'])
%
% now map the time averaged power input to the local shoreline coordinates:
xsl        = info.x_shoreline';
x_map      = x-xsl;
[xx,yy]    = meshgrid(x-mean(xsl),y);
if std(xsl)>=3
    fprintf('\tmapping to shoreline coordinates\n')
    Eavg_map   = griddata(x_map,yy,Eavg,xx,yy);
else
    Eavg_map   = Eavg;
end
%
alims = [-25 325 0 1500];
clims = [-1 1]*1e-2;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
xm = 2; ym = 2; pw = 10; ph = 3; ag=0.25; ps = [1.5*xm+pw, 2*ym+3*ag+2*ph];
ppos = [xm ym             pw ph];
cbpos = [xm ym+ph+ag 0.25*pw ag/1.5];
fig0  = figure('units','centimeters');
pos   = get(fig0,'position');
pos(3:4) = ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps])
pos0 = get(ax0,'position');
set(ax0,'position',pos0+[0 -0.75 0 0])
%
ax0  = axes('units','centimeters','position',ppos);
imagesc(ax0,y,x-mean(xsl),Eavg_map'),
% imagesc(ax0,y,x-mean(xsl),Eavg'), 
caxis(ax0,clims)
colormap(ax0,cm)
ylabel(ax0,'$y$ [m]','interpreter','latex')
xlabel(ax0,'$x$ [m]','interpreter','latex')
set(ax0,'tickdir','out','ticklabelinterpreter','latex','ylim',-25+[0 200],'ydir','normal','xdir','reverse','xticklabel',[])
%
% make colorbar
cx01 = axes('units','centimeters','position',cbpos);
imagesc(cx01,clrs,0,reshape(cm,1,256,3))
ylabel(cx01,{sprintf('$\left<u_\mathrm{rot}\cdot F_\mathrm{br}\right>_$ [m$^2$/s$^{3}$]'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',12,'yaxislocation','right')
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_power_input_time_averaged_vs_shoreline_coords.pdf'])
%
% compute k_y spectra on the time-averaged and shoreline mapped power input
dy = info.dy;
if ~mod(Ny,2)
    inyq = 1;% don't double the nyquist frequency
    stop = Ny/2+1;
else
    inyq = 0;%
    stop = (Ny+1)/2;
end
ky = [0:stop-1]'/(dy*Ny);
Eavg_map_spec = detrend(Eavg_map);
Eavg_map_spec = fft(Eavg_map_spec);
Eavg_map_spec = Eavg_map_spec(1:stop,:);
Eavg_map_spec = real(Eavg_map_spec.*conj(Eavg_map_spec));
Eavg_map_spec(2:end-inyq,:) = 2*Eavg_map_spec(2:end-inyq,:);
%
% average over 5-ky bins
Nf = 5;
ff = hamming(Nf)./hamming(Nf);
Eavg_map_spec= conv2(Eavg_map_spec,ff,'same');
% now cross-shore average into 10m bins (dof~2*20, less the cross-shore correlation)
db    = 25;
xbins = 0:db:150;
Nb    = length(xbins);
Eavg_map_spec_bin = nan(stop,Nb);
figure,
cm2 = cmocean('thermal',Nb);
for ii = 1:Nb
    inbin  =  find( (x-mean(xsl))>=(xbins(ii)-db/2) & (x-mean(xsl))<(xbins(ii)+db/2));
    Eavg_map_spec_bin(:,ii) = nanmean(Eavg_map_spec(:,inbin),2);
    loglog(ky,Eavg_map_spec_bin(:,ii),'color',cm2(ii,:)), hold on
end
% $$$ % what is the domain average of Prot vs. time?
% $$$ t              = rot_vel.t;
% $$$ Prot_SZavg = squeeze(nanmean(nanmean(Prot,1),2));
% $$$ % wave-average?
% $$$ % estimate ky-spectra of Prot?
