function info = estimate_FUNWAVE_velocity_helmholtz_decomposition(info);
%
% info = estimate_FUNWAVE_velocity_helmholtz_decomposition(info);
%
% calculate the helmholtz decompostion, u = u_rot + u_irr, where u_rot
% is the rotational component of the flow with stream function PSI,
% and u_irr is the irrotational component with potential function
% PHI. Analysis is limited to the region:
%       75 m <= x <=  525 m
%     -200 m <= y <= 1300 m

% get the input (x,y,h) grid file
load(info.bathyFile)
% switch from h=-z_bottom to h=depth
h = -h;
% preserve the original grid
x0=x;
y0=y;
h0=h;
%
% depth threshold
threshold = 0.01;
xsl       = info.x_shoreline;
avg_xsl   = round(mean(xsl));
mask0     = info.mask;
%
% full alongshore domain
iX = find(x0>=75  & x0<=525);
nx = length(iX);
ny = length(y0);
subDomain = [1 ny iX(1) iX(end)];
%
% [u,v,x,y,h,t,mask] = load_funwave_u_v(info.rootMat,info.rootName,info.bathyFile,subDomain);
% map bottom points to eta points
x = x(subDomain(3):subDomain(4));
y = y(subDomain(1):subDomain(2));
h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
             h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
nx = length(x);
ny = length(y);
dy = median(diff(y));
dx = median(diff(x));
% 
[xx,yy] = meshgrid(x,y);
%
% keep record of time and eta with velocity decomp
eta0  = [];
t0    = [];
% load (u,v) data and estimate (Urot,Vrot) for each file
files1 = dir([info.rootMat,info.rootName,'eta*.mat']);
files2 = dir([info.rootMat,info.rootName,'u*.mat']);
files3 = dir([info.rootMat,info.rootName,'v*.mat']);
Nf     = length(files1);
ROT    = matfile([info.rootMat,info.rootName,'rotational_velocity.mat'],'writable',true);
IRR    = matfile([info.rootMat,info.rootName,'irrotational_velocity.mat'],'writable',true);
ROT.x  = x;
ROT.y  = y;
ROT.h  = h;
%
% for averages
inds   = 0;
Uex    = 0;
U_avg  = 0;
V_avg  = 0;
Urot_avg=0;
Vrot_avg=0;
Erot    = 0;
%
%
% make a video of [vid1 = vorticity, vid2 = PSI]
alims = [75 400 0 1500];
clims0 = [-0.2 0.2];
clrs0  = clims0(1):diff(clims0)/255:clims0(2);
% $$$ clims1 = [-0.2 0.2];
% $$$ clrs1  = clims1(1):diff(clims1)/255:clims1(2);
% $$$ cm1    = cmocean('curl');
cm0    = cmocean('balance');
[fig0,ax0,ax00,cx0,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
% $$$ [fig1,ax1,ax00,cx1,ps,ppos,pos] = get_1panel_video_figure_info(alims);
% $$$ delete(ax00)
%
% generate video object
dat1 = matfile([files1(1).folder,'/',files1(1).name]);
t    = dat1.t;
t0   = t(1);
dt   = median(diff(t));
Navg   = 30;
Nplt   = 10;
vidName= sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'vorticity_avg%ds'],round(Navg*dt));
% $$$ vidName1= sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'rotational_streamfunction_avg%ds'],round(Navg*dt));
%
vid0 = VideoWriter(vidName,'Motion JPEG AVI');
vid0.Quality  = 100;
vid0.FrameRate= 5;
open(vid0)
%
% $$$ vid2 = VideoWriter(vidName2,'Motion JPEG AVI');
% $$$ vid2.Quality  = 100;
% $$$ vid2.FrameRate= 5;
% $$$ open(vid2)
%
for ii=1:Nf
    fprintf('loading u from: %s \n', files2(ii).name);
    if exist('subDomain','var')
        dat1 = matfile([files1(ii).folder,'/',files1(ii).name]);
        dat2 = matfile([files2(ii).folder,'/',files2(ii).name]);
        dat3 = matfile([files3(ii).folder,'/',files3(ii).name]);
        eta  = dat1.eta(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        t    = dat1.t;
        u    = dat2.u(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        v    = dat3.v(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
    else
        load([files1(ii).folder,'/',files1(ii).name]);        
        load([files2(ii).folder,'/',files2(ii).name],'u');
        load([files3(ii).folder,'/',files3(ii).name],'v');        
    end
% $$$     eta0 = cat(3,eta0,eta);        
% $$$     t0   = cat(1,t0,t);
    H    = h+eta;
    mask = (H>threshold);
    u(~mask) = nan;
    v(~mask) = nan;
    nt   = length(t);
    if ii==1
        dt = mean(diff(t));
    end
    II   = 0;
    for jj=1:nt;
        II = II+1;
        U  = u(:,:,jj);
        V  = v(:,:,jj);
        [omega,~] = curl(xx,yy,U,V);
        %
        [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(U,V,dx,dy);
        VORT(:,:,II)= omega;
        PSI(:,:,II) = psi;
        Urot(:,:,II)= u_psi;
        Vrot(:,:,II)= v_psi;
        PHI(:,:,II) = phi;
        Uirr(:,:,II)= u_phi;
        Virr(:,:,II)= v_phi;
        err(:,:,II) = sqrt( (U-(u_psi+u_phi)).^2 + (V-(v_psi+v_phi)).^2 );
        if II>Navg & mod(II,Nplt)==0
            i2plot = II-(Navg+1):II;
            %
            % MASKavg = min(mask(:,:,inds),[],3);
            VORTavg = nanmean(VORT(:,:,i2plot),3);% .*MASKavg;
            %
            figure(fig0)
            imagesc(ax0,y,x,VORTavg'), hold(ax0,'on'), plot(ax0,y,info.x_shoreline,'--k') 
            caxis(ax0,clims0)
            colormap(ax0,cm0)
            ylabel(ax0,'$y$ [m]','interpreter','latex')
            xlabel(ax0,'$x$ [m]','interpreter','latex')
            title(ax0,sprintf(' $t = %2.1f$ [min]~~~~', (mean(t(i2plot))-t0)/60),'interpreter','latex','horizontalalignment','right')
            set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
            %
            % make colorbar
            imagesc(cx0,clrs0,0,reshape(cm0,1,256,3))
            ylabel(cx0,{'$\left<\nabla_\mathrm{\scriptscriptstyle H} \times u_\mathrm{rot}\right>$ [s$^{-1}$]~~'},'interpreter','latex','rotation',0,'horizontalalignment','right')
            set(cx0,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
            %
            % get+write frame
            set(fig0,'position',pos);
            frame = getframe(fig0);
            if vid0.FrameCount==0
                fsize = size(frame.cdata,1,2);
            elseif any(size(frame.cdata,1,2)~=fsize)
                frame.cdata = imresize(frame.cdata,fsize);
            end
            writeVideo(vid0,frame)
        end
    end
    VORT(~mask) = nan;
    PSI (~mask) = nan;
    Urot(~mask) = nan;
    Vrot(~mask) = nan;
    PHI (~mask) = nan;
    Uirr(~mask) = nan;
    Virr(~mask) = nan;
    err(~mask)  = nan;
    %
    inds = inds(end)+[1:nt];
    ROT.t(inds,1)            = t;
    ROT.eta(1:ny,1:nx,inds)  = eta;
    ROT.mask(1:ny,1:nx,inds) = mask;
    ROT.VORT(1:ny,1:nx,inds) = VORT;
    ROT.PSI (1:ny,1:nx,inds) = PSI;
    ROT.Urot(1:ny,1:nx,inds) = Urot;
    ROT.Vrot(1:ny,1:nx,inds) = Vrot;
    ROT.err (1:ny,1:nx,inds) = err;
    ROT.t(inds,1)            = t;    
    IRR.eta(1:ny,1:nx,inds)  = eta;
    IRR.mask(1:ny,1:nx,inds) = mask;
    IRR.PHI (1:ny,1:nx,inds) = PHI;
    IRR.Uirr(1:ny,1:nx,inds) = Uirr;
    IRR.Virr(1:ny,1:nx,inds) = Virr;
    IRR.err (1:ny,1:nx,inds) = err;    
    %
    % need to map these to a shoreline following x-coord before averagting
    absU_avg_t = nanmean(abs(Urot),3);
    avg_absU = 0*x;
    for kk = 1:ny
        tmp1 = interp1( (x-xsl(kk)), absU_avg_t(kk,:), x-avg_xsl);
        avg_absU = avg_absU+tmp1;
    end
    Uex = Uex + avg_absU/ny;% nansum(nansum(abs(Urot).*H*dt*dy,1),3)./nansum(nansum(H*dy*dt,1),3)
    U_avg    = U_avg + nansum(u,3)/nt;
    V_avg    = V_avg + nansum(v,3)/nt;        
    Urot_avg = Urot_avg + nansum(Urot,3)/nt;
    Vrot_avg = Vrot_avg + nansum(Vrot,3)/nt;
    Erot  = Erot + nansum(Urot.^2 + Vrot.^2,3)/nt;
    clear VORT PSI Urot Vrot PHI Uirr Virr omega err
end
% quick stats
Uex = Uex/Nf;
U_avg = U_avg/Nf;
V_avg = V_avg/Nf;
Urot_avg = Urot_avg/Nf;
Vrot_avg = Vrot_avg/Nf;
% calculate exchange due to mean flow
Uex_avg  = 0*x;
for kk=1:ny
    tmp1 = interp1( (x-xsl(kk)), abs(Urot_avg(kk,:)), x-avg_xsl);
    Uex_avg = Uex_avg + tmp1;
end
Uex_avg  = Uex_avg/ny;
% calculate rotational energy 
Erot = Erot/Nf;
ROT.Uex  = Uex;
ROT.Uex_avg=Uex_avg;
ROT.Urot_avg = Urot_avg;
ROT.Vrot_avg = Vrot_avg;
ROT.Erot  = Erot;
Eeddy    = Erot - Urot_avg.^2 - Vrot_avg.^2;
ROT.Eeddy= Eeddy;
%
% make simple stats plot
fig0 = figure;
p0   = plot(x,Uex,'-k',x,Uex_avg,'-b',x,Uex-Uex_avg,'-r','linewidth',2);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$U_\mathrm{\scriptscriptstyle{EX}}$ [m/s]','interpreter','latex')
set(gca,'xlim',[75 525],'ticklabelinterpreter','latex','tickdir','out')
f0l1 = legend([p0(1),p0(2),p0(3)]','$U_\mathrm{\scriptscriptstyle EX}$','$\left<u\right>_\mathrm{\scriptscriptstyle EX}$','$u''_\mathrm{\scriptscriptstyle EX}$');
set(f0l1,'location','northeast','interpreter','latex')
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ',[info.rootSim,filesep,'figures']])
end
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'Uex_exchange_velocity.pdf'])
%
vel_range = max(range(Urot_avg(:)),range(Vrot_avg(:)));
clims    = vel_range/2*[-1 1];
clrs     = clims(1):diff(clims)/255:clims(2);
cm       = cmocean('balance');
fig1 = figure;
a1 = subplot(1,2,1)
imagesc(x,y,Urot_avg,clims),colormap(cm),caxis(clims)
xlabel(a1,'crosshore [m]','interpreter','latex')
ylabel(a1,'alongshore [m]','interpreter','latex')
title('$\left<U_\mathrm{rot}\right>$','interpreter','latex')
set(a1,'xlim',[75 500],'ticklabelinterpreter','latex','tickdir','out','ydir','normal')
a2 = subplot(1,2,2);
imagesc(x,y,Vrot_avg,clims),colormap(cm),caxis(clims)
xlabel(a2,'crosshore [m]','interpreter','latex')
ylabel(a2,'alongshore [m]','interpreter','latex')
title(a2,'$\left<V_\mathrm{rot}\right>$','interpreter','latex')
set(a2,'xlim',[75 500],'ticklabelinterpreter','latex','tickdir','out','ydir','normal')
a3 = axes('position',[0.85 0.15 0.025 0.2]);
imagesc(0,clrs,reshape(cm,256,1,3))
ylabel(a3,'[m/s]~~','rotation',0,'interpreter','latex','horizontalalignment','right')
set(a3,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xtick',[])
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'Urot_Vrot_time_averaged.pdf'])
%
%
vel_range = max(range(U_avg(:)),range(V_avg(:)));
clims    = vel_range/2*[-1 1];
clrs     = clims(1):diff(clims)/255:clims(2);
cm       = cmocean('balance');
fig1v2 = figure;
a1 = subplot(1,2,1)
imagesc(x,y,U_avg,clims),colormap(cm),caxis(clims)
xlabel(a1,'crosshore [m]','interpreter','latex')
ylabel(a1,'alongshore [m]','interpreter','latex')
title('$\left<U\right>$','interpreter','latex')
set(a1,'xlim',[75 500],'ticklabelinterpreter','latex','tickdir','out','ydir','normal')
a2 = subplot(1,2,2);
imagesc(x,y,V_avg,clims),colormap(cm),caxis(clims)
xlabel(a2,'crosshore [m]','interpreter','latex')
ylabel(a2,'alongshore [m]','interpreter','latex')
title(a2,'$\left<V\right>$','interpreter','latex')
set(a2,'xlim',[75 500],'ticklabelinterpreter','latex','tickdir','out','ydir','normal')
a3 = axes('position',[0.85 0.15 0.025 0.2]);
imagesc(0,clrs,reshape(cm,256,1,3))
ylabel(a3,'[m/s]~~','rotation',0,'interpreter','latex','horizontalalignment','right')
set(a3,'ticklabelinterpreter','latex','fontsize',10,'tickdir','out','xtick',[])
exportgraphics(fig1v2,[info.rootSim,filesep,'figures',filesep,info.rootName,'U_V_time_averaged.pdf'])
%
%
Ly   = info.Ly;
Et   = nansum(Erot*dy,1)/Ly;
Em   = nansum(Urot_avg.^2 + Vrot_avg.^2,1)/Ly;
Ek   = nansum(Eeddy*dy,1)/Ly;
fig2 = figure;
p0   = plot(x,sqrt(Et),'-k',x,sqrt(Em),'-b',x,sqrt(Ek),'linewidth',2);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$\mathcal{E}^{1/2}$ [m/s]','interpreter','latex')
f0l1 = legend([p0(1),p0(2),p0(3)]','$u$','$\left<u\right>$','$u''$');
set(f0l1,'location','northeast','interpreter','latex')
set(gca,'xlim',[75 525],'ticklabelinterpreter','latex','tickdir','out')
exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'Erot_eddy_kinetic_energy.pdf'])
%
%
close all
%
info.rotVelFile = ROT.Properties.Source;% [info.rootMat,info.rootName,'rotational_velocity.mat'];
info.irrVelFile = IRR.Properties.Source;% [info.rootMat,info.rootName,'irrotational_velocity.mat'];
save(info.fileName,'-struct','info')
