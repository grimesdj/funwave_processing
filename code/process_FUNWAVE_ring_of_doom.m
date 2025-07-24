function info = process_FUNWAVE_ring_of_doom(info)

stationsFile = [info.rootMat,info.rootName,'gauges.mat'];
load(stationsFile)
% estimate PUV wave stats at each location
%
% get the instruments around the ROD, load locations
load('/home/derek/projects/ShortCrests/ROD/mat_data/FRF_instrument_locations_2013.mat','ADVs','rodX','rodY')
[rodXrot,rodYrot] = rotate_FRF_coords(rodX,rodY,info.angle);
[advX   ,   advY] = rotate_FRF_coords(ADVs(:,1),ADVs(:,2),info.angle);
%
% find instruments within 5m radius of (x,y)_ROD
D      = sqrt( ([stations.x]-rodXrot).^2 + ([stations.y]-rodYrot).^2 );
iROD   = find(D<=5);
[~,i0] = min(D);
%
% now get the coordinates of each ADV on the ROD.
theta0   = 0;
dtheta   = 360/14;
thetaADV = theta0 + [0:13]*dtheta;% :dtheta:360-dtheta;
rROD     = 2.5;
% get coordinates for ADVs at r=2.5m from ring center
xROD_ADV = stations(i0).x+rROD*cosd(thetaADV);
yROD_ADV = stations(i0).y+rROD*sind(thetaADV);
%
% for each ROD_ADV, either use one of the stations or a cluster to interpolate
ROD_MOD  = struct([]);
ROD_MOD(1).U = stations(i0).u;%interp1(stations(i0).t,stations(i0).u,time);
ROD_MOD(1).V = stations(i0).v;%interp1(stations(i0).t,stations(i0).v,time);
ROD_MOD(1).H = stations(i0).h+stations(i0).eta;%interp1(stations(i0).t,stations(i0).h+stations(i0).eta,time);
time      = stations(i0).t;
i1       = iROD(iROD~=i0);
DD       = sqrt( ([stations(i1).x]'-xROD_ADV).^2 + ([stations(i1).y]'-yROD_ADV).^2 );
L        = 2*rROD*tan(dtheta/2*pi/180);
A        = 0.5*L*rROD*14;
for jj=1:14
    distances   = DD(:,jj);
    if any(distances==0)
        coloc = find(distances==0);
        ROD_MOD(1).x(jj) = stations(i1(coloc)).x;
        ROD_MOD(1).y(jj) = stations(i1(coloc)).y;
        ROD_MOD(1).h(jj) = stations(i1(coloc)).h;
        ROD_MOD(1).t     = time;        
        ROD_MOD(1).eta(:,jj) = stations(i1(coloc)).eta;
        ROD_MOD(1).u  (:,jj) = stations(i1(coloc)).u;
        ROD_MOD(1).v  (:,jj) = stations(i1(coloc)).v;
    else
        %
        instruments = find(distances<=sqrt(info.dx^2+info.dy^2));
        interpWeight= distances(instruments).^(-2)./sum(distances(instruments).^(-2));
        ROD_MOD.x(jj) = [stations(i1(instruments)).x]*interpWeight;
        ROD_MOD.y(jj) = [stations(i1(instruments)).y]*interpWeight;
        ROD_MOD.h(jj) = [stations(i1(instruments)).h]*interpWeight;                
        ROD_MOD.t     = time;        
        ROD_MOD.eta(:,jj) = interp1([stations(i1(instruments)).t  ]*interpWeight,...
                                    [stations(i1(instruments)).eta]*interpWeight,time);
        ROD_MOD.u  (:,jj) = interp1([stations(i1(instruments)).t  ]*interpWeight,...
                                    [stations(i1(instruments)).u  ]*interpWeight,time);
        ROD_MOD.v  (:,jj) = interp1([stations(i1(instruments)).t  ]*interpWeight,...
                                    [stations(i1(instruments)).v  ]*interpWeight,time);
    end
        ROD_MOD.Ur (:,jj) = ROD_MOD.u(:,jj)*cosd(thetaADV(jj)) + ROD_MOD.v(:,jj)*sind(thetaADV(jj));
        ROD_MOD.Ua (:,jj) =-ROD_MOD.u(:,jj)*sind(thetaADV(jj)) + ROD_MOD.v(:,jj)*cosd(thetaADV(jj));
        ROD_MOD.ua (:,jj) =-ROD_MOD.Ua(:,jj)*sind(thetaADV(jj));% x-component of azimuthal vel
        ROD_MOD.va (:,jj) = ROD_MOD.Ua(:,jj)*cosd(thetaADV(jj));% y-component of azimuthal vel        
        ROD_MOD.ur (:,jj) = ROD_MOD.Ur(:,jj)*cosd(thetaADV(jj));% ... radial
        ROD_MOD.vr (:,jj) = ROD_MOD.Ur(:,jj)*sind(thetaADV(jj));% ... radial        
end
ROD_MOD.vort = sum(L/A*ROD_MOD.Ua.*(ROD_MOD.h+ROD_MOD.eta),2)./sum(ROD_MOD.h+ROD_MOD.eta,2);
ROD_MOD.div  = sum(L/A*ROD_MOD.Ur.*(ROD_MOD.h+ROD_MOD.eta),2)./sum(ROD_MOD.h+ROD_MOD.eta,2);
%
%
% load the ROD data structure
rodFile = sprintf('/home/derek/projects/ShortCrests/ROD/data/RD_2013%s.mat',info.dateTime(1:4));
load(rodFile)
hr     = str2num(info.dateTime(5:6));
%
% load the local vorticity from funwave output
load(info.rotVelFile,'VORT','t')
% load grided (x,y) locations to find index nearest the ROD
load(info.bathyFile,'x','y','h')
[~,iX] = min(abs(x-rodXrot));
[~,iY] = min(abs(y-rodYrot));
indROD = sub2ind(size(h),iY,iX);
depROD =-h(indROD);
ROD_MOD.time_local = t;
ROD_MOD.vort_local = squeeze(VORT(iY,iX,:));
ROD_MOD.dep_local  = depROD;
%
%
%
ll = 100;
fig0 = figure;
plot(ROD_MOD.x,ROD_MOD.y,'xk',...
     ROD_MOD.x + [0;1]*ROD_MOD.ua(ll,:), ROD_MOD.y + [0;1]*ROD_MOD.va(ll,:),'-r',...
     ROD_MOD.x + ROD_MOD.ua(ll,:) + [0;1]*ROD_MOD.ur(ll,:), ROD_MOD.y + ROD_MOD.va(ll,:) + [0;1]*ROD_MOD.vr(ll,:),'-c',...
     ROD_MOD.x + [0;1]*ROD_MOD.u (ll,:), ROD_MOD.y + [0;1]*ROD_MOD.v (ll,:),':k')
axis equal
grid on
xlabel('crosshore [m]','interpreter','latex')
ylabel('alongshore [m]','interpreter','latex')
text(mean(ROD_MOD.x)-1,mean(ROD_MOD.y),{ sprintf('$\\bar{\\omega} =$ %1.4f~s$^{-1}$',ROD_MOD.vort(ll));...
                                         sprintf('$\\bar{u} =$ %1.4f~m/s',mean(ROD_MOD.u(ll,:)));...
                                         sprintf('$\\bar{v} =$ %1.4f~m/s',mean(ROD_MOD.v(ll,:)))},'interpreter','latex','fontsize',12)
set(gca,'ticklabelinterpreter','latex','tickdir','out','xlim',mean(ROD_MOD.x)+[-4 4],'ylim',mean(ROD_MOD.y)+[-4 4],'fontsize',14)
pause(0.1)
%
if ll
    if ~exist([info.rootSim,filesep,'figures'],'dir')
        eval(['!mkdir ',[info.rootSim,filesep,'figures']])
    end
    exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_velocity_components_example.pdf'])
    %
end
%
fig1 = figure;
VORTfactor = (std(RD(hr).vort,'omitnan')./std(ROD_MOD.vort,'omitnan'))
VORTlocalFactor = (std(RD(hr).vort,'omitnan')./std(ROD_MOD.vort_local,'omitnan'))
plot(RD(hr).t_sec-RD(hr).t_sec(1),RD(hr).vort,'.k',ROD_MOD.t-ROD_MOD.t(1),ROD_MOD.vort*VORTfactor,'.r',ROD_MOD.time_local-ROD_MOD.time_local(1),ROD_MOD.vort_local*VORTlocalFactor,'.b')
xlabel('time [s]','interpreter','latex')
ylabel('$\bar{\omega}$ [s$^{-1}$]','interpreter','latex')
text(50,0.6,sprintf('$%f\\times\\bar\\omega_\\mathrm{mod}$',VORTfactor),'color','r')
text(50,0.5,sprintf('$%f\\times\\omega_\\mathrm{mod}$',VORTlocalFactor),'color','b')
set(gca,'ticklabelinterpreter','latex','tickdir','out','xlim',[0 3e3],'fontsize',14)
exportgraphics(fig1,sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_vorticity_obs_mod.pdf'],VORTfactor))
%
%
%
[Sww      ,f      ]=welch_method(ROD_MOD.vort(1:end),0.25,25,0.5);
[Sww_local,f_local]=welch_method(ROD_MOD.vort_local(1:end),1,12,0.5);
fig2 = figure;
semilogx(RD(hr).fm,RD(hr).Svort,'-k',f,Sww*VORTfactor^2,'-r',f_local,Sww_local*VORTlocalFactor^2,'-b')
xlabel('frequency [Hz]','interpreter','latex')
ylabel('$S_{\omega}$ [s$^{-2}$/Hz]','interpreter','latex')
text(1/10,0.15,sprintf('$%f\\times S_{\\bar\\omega_\\mathrm{mod}}$',VORTfactor^2),'color','r')
text(1/10,0.1,sprintf('$%f\\times S_{\\omega_\\mathrm{mod}}$',VORTlocalFactor^2),'color','b')
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',14)
exportgraphics(fig2,sprintf([info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_vorticity_spectra_obs_mod.pdf'],VORTfactor^2))
%
ROD_MOD.VORTfactor = VORTfactor;
ROD_MOD.freq = f;
ROD_MOD.Svort= Sww;
ROD_MOD.VORTlocalFactor = VORTlocalFactor;
ROD_MOD.freq_local = f_local;
ROD_MOD.Svort_local= Sww_local;
%
%
modRODfile = [info.rootMat,info.rootName,'ROD_statistics.mat'];
info.RODfile = modRODfile;
%
save(info.fileName,'-struct','info')
ROD_OBS = RD(hr);
save(modRODfile,'ROD_MOD','ROD_OBS')