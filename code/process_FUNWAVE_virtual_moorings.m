function info = process_FUNWAVE_virtual_moorings(info)
%
% USAGE: info = process_FUNWAVE_virtual_moorings(info)
%
% load virtual mooring data and estimate wave/current/vorticity stats
% !!!need to align the ADVs with the modeled moorings!!!

g       = 9.8;
% get the station locations
if isfield(info,'stationFile')
    stationID = fopen(info.stationFile);
elseif isfield(info,'gaugeFile')
    stationID = fopen(info.gaugeFile);
else
    fprintf('\n\t!!!!! No Gauge/Station File: %s !!!!!\n',info.runName)
    return
end
stationXY = textscan(stationID,'%f %f','delimiter','\t');
iX        = stationXY{1};
iY        = stationXY{2};
Ns        = size(iX,1);
Is        = [1:Ns]';

% load grided (x,y) locations
load(info.bathyFile,'x','y','h')
Xs = x(iX);
Ys = y(iY);
inds = sub2ind(size(h),iY,iX);
dep  =-h(inds);

% download the station files
stationFiles = dir([info.rootSim,filesep,'output',filesep,'sta_*']);
Nf           = length(stationFiles);
if Nf~=Ns
    disp('different number of stations (x,y) and station files')
    return
end
%
stations = struct([]);
for ii = 1:Nf
    fid = fopen([stationFiles(ii).folder,filesep,stationFiles(ii).name]);
    data= textscan(fid,'%f %f %f %f','delimiter','\t');
    tmp = split(stationFiles(ii).name,'_');
    ID  = str2num(tmp{2});
    ind = find(ID==Is);

    % data is not sampled at regular intervals
    t = data{1};
    t = t-(t(end)-3600);% only look at 1-hour
    [tu,uni] = unique(t);
    if ii==1
        dt   = 0.125;
        fprintf('\n\t!!!warning: assuming dt = %f s !!!!\n',dt)
        time = [0:dt:t(end)-dt]';
        Nt   = length(time);
        % stuff for computing spectra
        fs = 1./dt;
        %     number of samples in average
        Na = 2^floor(log2(Nt));%2^12;
        %     number of samples in each ensemble
        Ne = 2^11;
        olap = 1/2;
        chnks = (Na-Ne*olap-1)/(Ne*(1-olap));
    end
    % extract fields
    eta= interp1(tu,data{2}(uni),time);
    u  = interp1(tu,data{3}(uni),time);
    v  = interp1(tu,data{4}(uni),time);
    x  = Xs(ind);
    y  = Ys(ind);
    h  = dep(ind);
    % compute spectra
    [Suu,fq]=welch_method(u(1:Na),dt,chnks,0.5);
    [Svv,fq]=welch_method(v(1:Na),dt,chnks,0.5);
    [See,fq]=welch_method(eta(1:Na),dt,chnks,0.5);
    [Suv,fq]=welch_cospec(u(1:Na),v(1:Na),dt,chnks,0.5);        
    [Spu,fq]=welch_cospec(eta(1:Na),u(1:Na),dt,chnks,0.5);
    [Spv,fq]=welch_cospec(eta(1:Na),v(1:Na),dt,chnks,0.5);
    % estimate bulk wave statistics
    fq = fq(2:end);
    Suu= Suu(2:end,:);
    Svv= Svv(2:end,:);
    Suv= Suv(2:end,:);
    See= See(2:end,:);
    Spu= Spu(2:end,:);
    Spv= Spv(2:end,:);
    %
    % convert (u,v) to surface elevation spectra
    alpha  = 0.531;
    om     = 2*pi*fq;
    lom=length(om);
    g      = 9.81;
    k      = wavenumber(om,h);
    cU2eta = (om./(g*k)).*cosh(k.*h)./cosh(k.*h.*(1-alpha));
    % kill frequencies above 1Hz?
    Beta   = 0.5*(1-tanh( (fq-0.5)/0.05 ));
    % surface transform spectra
    SeUU = Suu.*(cU2eta.*Beta).^2;
    SeVV = Svv.*(cU2eta.*Beta).^2;
    SeUV = Suv.*(cU2eta.*Beta).^2;
    SePP = See;
    SePU = Spu.*cU2eta.*Beta;
    SePV = Spv.*cU2eta.*Beta;
    % 1c) calculate kuik parameters
    coSue =   SePU;
    coSve =   SePV;
    coUVe =   SeUV;
    r2d = 180/pi;
    %
    % get (a1,b1) & (a2,b2)
    a1      = coSue ./ sqrt( SePP .* ( SeUU + SeVV ) );
    b1      = coSve ./ sqrt( SePP .* ( SeUU + SeVV ) );
    dir1    = r2d* ( atan2(b1,a1) );          
    %
    % average over wind-wave band
    df= fq(2)-fq(1);
    I = find(fq>=1/20 & fq<=1/4);
    [~,Ip] = max(SePP,[],1);
    m0 = nansum(SePP(I,:)*df);
    m1 = nansum(fq(I).*SePP(I,:)*df);        
    ma1= nansum(a1(I,:).*SePP(I,:)*df,1)./m0;
    mb1= nansum(b1(I,:).*SePP(I,:)*df,1)./m0;
    mdir1=r2d*atan2(mb1,ma1);
    mspread1 = r2d*sqrt(abs(2*(1-(ma1.*cos(mdir1/r2d) + mb1.*sin(mdir1/r2d)))));
    %
    % get a2 b2
    a2 = (SePP - SeVV) ./ (SePP + SeVV);
    b2 = 2 .* coUVe ./ ( SePP + SeVV );
    %
    ma2= nansum(a2(I,:).*SePP(I,:)*df,1)./m0;
    mb2= nansum(b2(I,:).*SePP(I,:)*df,1)./m0;
    mdir2=r2d/2*atan2(mb2,ma2);
    mspread2 = r2d*sqrt(abs(0.5-0.5*(ma2.*cos(2.*mdir1/r2d)+mb2.*sin(2.*mdir1/r2d))));
    %
    Hs = 4*sqrt(m0);
    Tm = m0./m1;
    % estimate rotational/irrotational current magnitude
    Stot_mod = Suu + Svv;
    Sirr_mod = min(See*(g/h),Stot_mod);
    Srot_mod = Stot_mod - Sirr_mod;
    %
    % get the vlf/ig band
    Imod = find(fq>0.004 & fq<0.03);
    Urot_mod = sqrt(sum(Srot_mod(Imod)*df));
    Rrot_mod = Stot_mod/Sirr_mod;
    %
    % log stats
    IDs(ii)= ID;
    stations(ii).t  = time;
    stations(ii).Nt = Nt;
    stations(ii).eta= eta;
    stations(ii).u  = u;
    stations(ii).v  = v
    stations(ii).x  = x;
    stations(ii).y  = y;
    stations(ii).h  = h;
    stations(ii).freq= fq;    
    stations(ii).See= See;
    stations(ii).Suu= Suu;
    stations(ii).Svv= Svv;
    stations(ii).Hs = Hs;
    stations(ii).Tm = Tm;
    stations(ii).a1b1= [ma1, mb1];
    stations(ii).a2b2= [ma2, mb2];
    stations(ii).dire= [mdir1,mdir2];
    stations(ii).spread=[mspread1,mspread2];
    stations(ii).Sirr = Sirr_mod;
    stations(ii).Srot = Srot_mod;    
    stations(ii).Urot = Urot_mod;
    stations(ii).Rrot = Rrot_mod;
end
stations(IDs)=stations;
%
stationsFile = [info.rootMat,info.rootName,'gauges.mat'];
save(stationsFile,'stations');
info.virtualMooringFile = stationsFile;
save(info.fileName,'-struct','info')
%
%
% 1) pair the stations with ADVs and
% 2) plot some statistics of velocity and sea-surface fluctuations


return
%
% i'm moving the below code to a new file: process_FUNWAVE_ring_of_doom.m
%
% $$$ if any(Nt~=min(Nt))
% $$$     disp('different length time vectors for stations')
% $$$     keyboard
% $$$ end
%
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
ROD_MOD(1).U = interp1(stations(i0).t,stations(i0).u,time);
ROD_MOD(1).V = interp1(stations(i0).t,stations(i0).v,time);
ROD_MOD(1).H = interp1(stations(i0).t,stations(i0).h+stations(i0).eta,time);
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
        ROD_MOD.ua (:,jj) =-ROD_MOD.Ua(:,jj)*sind(thetaADV(jj));
        ROD_MOD.va (:,jj) = ROD_MOD.Ua(:,jj)*cosd(thetaADV(jj));        
        ROD_MOD.ur (:,jj) = ROD_MOD.Ur(:,jj)*cosd(thetaADV(jj));
        ROD_MOD.vr (:,jj) = ROD_MOD.Ur(:,jj)*sind(thetaADV(jj));        
end
ROD_MOD.vort = sum(L/A*ROD_MOD.Ua.*(ROD_MOD.h+ROD_MOD.eta),2)./sum(ROD_MOD.h+ROD_MOD.eta,2);
ROD_MOD.div  = sum(L/A*ROD_MOD.Ur.*(ROD_MOD.h+ROD_MOD.eta),2)./sum(ROD_MOD.h+ROD_MOD.eta,2);
%
%
fig0 = figure;
for ll=1% :100;
plot(ROD_MOD.x,ROD_MOD.y,'xk',...
     ROD_MOD.x + [0;1]*ROD_MOD.ua(2e3+ll,:), ROD_MOD.y + [0;1]*ROD_MOD.va(2e3+ll,:),'-r',...
     ROD_MOD.x + ROD_MOD.ua(2e3+ll,:) + [0;1]*ROD_MOD.ur(2e3+ll,:), ROD_MOD.y + ROD_MOD.va(2e3+ll,:) + [0;1]*ROD_MOD.vr(2e3+ll,:),'-c',...
     ROD_MOD.x + [0;1]*ROD_MOD.u (2e3+ll,:), ROD_MOD.y + [0;1]*ROD_MOD.v (2e3+ll,:),':k')
axis equal
grid on
xlabel('crosshore [m]','interpreter','latex')
ylabel('alongshore [m]','interpreter','latex')
text(mean(ROD_MOD.x)-1,mean(ROD_MOD.y),{ sprintf('$\\bar{\\omega} =$ %1.4f~s$^{-1}$',ROD_MOD.vort(2e3+ll));...
                                         sprintf('$\\bar{u} =$ %1.4f~m/s',mean(ROD_MOD.u(2e3+ll,:)));...
                                         sprintf('$\\bar{v} =$ %1.4f~m/s',mean(ROD_MOD.v(2e3+ll,:)))},'interpreter','latex','fontsize',12)
set(gca,'ticklabelinterpreter','latex','tickdir','out','xlim',mean(ROD_MOD.x)+[-4 4],'ylim',mean(ROD_MOD.y)+[-4 4],'fontsize',14)
pause(0.1)
%
if ll==1
    if ~exist([info.rootSim,filesep,'figures'],'dir')
        eval(['!mkdir ',[info.rootSim,filesep,'figures']])
    end
    exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_velocity_components_example.pdf'])
    %
end
end
%
%
% 
%
%
rodFile = sprintf('/home/derek/projects/ShortCrests/ROD/data/RD_2013%s.mat',info.dateTime(1:4));
load(rodFile)
hr     = str2num(info.dateTime(5:6));
fig1 = figure;
plot(RD(hr).t_sec-RD(hr).t_sec(1),RD(hr).vort,'.k',ROD_MOD.t(2e3:end)-ROD_MOD.t(2e3),ROD_MOD.vort(2e3:end)*15,'.r')
xlabel('time [s]','interpreter','latex')
ylabel('$\bar{\omega}$ [s$^{-1}$]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','xlim',[0 3e3],'fontsize',14)
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_vorticity_obs_modx15.pdf'])
%
[Sww,f]=welch_method(ROD_MOD.vort(2e3:end),dt,25,0.5);
fig2 = figure;
semilogx(RD(hr).fm,RD(hr).Svort,'-k',f,Sww*1e3,'-r')
xlabel('frequency [Hz]','interpreter','latex')
ylabel('$S_{\omega}$ [s$^{-2}$/Hz]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',14)
exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_vorticity_spectra_obs_modx1000.pdf'])
%
[Sdd,f]=welch_method(ROD_MOD.div (2e3:end),dt,25,0.5);
fig3 = figure; semilogx(RD(hr).fm,RD(hr).Sdiv,'-k',f,Sdd*1e5,'-r')
xlabel('frequency [Hz]','interpreter','latex')
ylabel('$S_\mathrm{div}$ [s$^{-2}$/Hz]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',14)
exportgraphics(fig3,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_divergence_spectra_obs_modx1000.pdf'])
%
%
% $$$ u = cosd(info.angle)*RD(hr).u + sind(info.angle)*RD(hr).v;
% $$$ v =-sind(info.angle)*RD(hr).u + cosd(info.angle)*RD(hr).v;
% $$$ %
% $$$ bins = [-3:0.05:3];
% $$$ [Hmod] = hist(ROD_MOD.u(1800:end,1),bins);
% $$$ [Hrod] = hist(u,bins);%RD(hr).u,bins);
% $$$ Hmod = Hmod./numel(ROD_MOD.u(1800:end,1));
% $$$ Hrod = Hrod./numel(RD(hr).u);
% $$$ figure,
% $$$ bar(bins,Hrod,'r')
% $$$ hold on,
% $$$ bar(bins,Hmod,'b','facealpha',0.25)
% $$$ %
% $$$ [Hmod] = hist(ROD_MOD.v(1800:end,1),bins);
% $$$ [Hrod] = hist(v,bins);%RD(hr).v,bins);
% $$$ Hmod = Hmod./numel(ROD_MOD.v(1800:end,1));
% $$$ Hrod = Hrod./numel(RD(hr).v);
% $$$ figure,
% $$$ bar(bins,Hrod,'r')
% $$$ hold on,
% $$$ bar(bins,Hmod,'b','facealpha',0.25)
%
[Suu,f]=welch_method(ROD_MOD.U(2e3:end,1),dt,25,0.5);
[Svv,f]=welch_method(ROD_MOD.V(2e3:end,1),dt,25,0.5);
[See,f]=welch_method(ROD_MOD.H(2e3:end,1),dt,25,0.5);
ROD_MOD.f   = f;
ROD_MOD.Suu = Suu;
ROD_MOD.Svv = Svv;
ROD_MOD.See = See;
%
% estimate rotational velocity spectra using Lippmann et al., 1999
g       = 9.8;
dep_obs = mean(RD(hr).p);
dep_mod = mean(ROD_MOD.H);
%
Stot_mod = Suu + Svv;
Sirr_mod = See*(g/dep_mod);
Srot_mod = Stot_mod - Sirr_mod;
%
Stot_obs = RD(hr).SUU + RD(hr).SVV;
Sirr_obs = RD(hr).SEee*(g/dep_obs);
Srot_obs = Stot_obs - Sirr_obs; Sirr_obs(Srot_obs<0)=Stot_obs(Srot_obs<0); Srot_obs(Srot_obs<0)=0;
%
% make some plots for Duke Marine Lab talk:
xm = 2;%cm
ym = 2;
ph = 5;
pw = 8;
ppos = [xm ym pw ph];
ps   = [1.25*xm+pw 1.25*ym+ph];
%
fig0 = figure('units','centimeters','color','w');
pos  = get(fig0,'position');
pos(3:4) = ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'inverthardcopy','off')
ax0  = axes('units','centimeters','position',ppos);
p0   = semilogx(RD(hr).fm, Stot_obs, '-k','linewidth',2);
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('[(m/s)$^2$/Hz]','interpreter','latex')
h0   = legend('$S_\mathrm{tot} = S_{u,u}+S_{v,v}$');
set(h0,'interpreter','latex','box','off','color','none','location','northwest')
set(ax0,'position',ppos,'xlim',[RD(hr).fm(2) 1],'ticklabelinterpreter','latex','xtick',[1e-3,1e-2,1e-1,1e0])
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_velocity_spectra.pdf'])
%
hold on,
p1 = semilogx(RD(hr).fm, Srot_obs,'-b',RD(hr).fm, Sirr_obs ,'-r','linewidth',2);
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('[(m/s)$^2$/Hz]','interpreter','latex')
h1   = legend([p0;p1],'$S_\mathrm{tot} = S_{u,u}+S_{v,v}$','$S_\mathrm{rot}$','$S_\mathrm{irr}$');
% $$$ set(h1,'interpreter','latex','box','off','color','none','location','northwest')
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'inverthardcopy','off')
set(ax0,'position',ppos,'xlim',[RD(hr).fm(2) 1],'ticklabelinterpreter','latex','xtick',[1e-3,1e-2,1e-1,1e0])
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_velocity_spectra_decomposition.pdf'])
%
%
fig4 = figure;
semilogx(RD(hr).fm,RD(hr).SUU+RD(hr).SVV,'-k',f,Suu+Svv,'-r',RD(hr).fm,Srot_obs,'--k',f,Srot_mod,'--r')
xlabel('frequency [Hz]','interpreter','latex')
ylabel('[(m/s)$^2$/Hz]','interpreter','latex')
%
Imod = find(f>0.004 & f<0.03);
Iobs = find(RD(hr).fm>0.004 & RD(hr).fm<0.03);
Eig_mod  =  sum(See(Imod))*(f(2)-f(1));
Uig_mod = (sum(Suu(Imod)) + sum(Svv(Imod)))*(f(2)-f(1));
Eig_obs  =  sum(RD(hr).SEee(Iobs))*(RD(hr).fm(2)-RD(hr).fm(1));
Uig_obs = (sum(RD(hr).SUU(Iobs)) + sum(RD(hr).SVV(Iobs)))*(RD(hr).fm(2)-RD(hr).fm(1));
%
Imod = find(f>0.04 & f<0.25);
Iobs = find(RD(hr).fm>0.04 & RD(hr).fm<0.25);
Hs_mod = 4*sqrt(sum(See(Imod).*(f(2)-f(1))));
Hs_obs = 4*sqrt(sum(RD(hr).SEee(Iobs).*(RD(hr).fm(2)-RD(hr).fm(1))));
%
Robs = (Uig_obs./Eig_obs)./(g./dep_obs);
Rmod = (Uig_mod./Eig_mod)./(g./dep_mod);
Urot_ig_mod = sqrt( Uig_mod .* (1-1/Rmod) );
Urot_ig_obs = sqrt( Uig_obs .* (1-1/Robs) );
text(0.3,4,{ sprintf('Mod:~~$H_\\mathrm{s} =$ %1.2f~m/s',Hs_mod);...
             sprintf('~~~~~~~~~$U_{\\psi,\\mathrm{\\scriptscriptstyle{IG}}} =$ %1.2f~m/s',Urot_ig_mod);...
             sprintf('Obs:~~$H_\\mathrm{s} =$ %1.2f~m/s',Hs_obs)
             sprintf('~~~~~~~~~$U_{\\psi,\\mathrm{\\scriptscriptstyle{IG}}} =$ %1.2f~m/s',Urot_ig_obs)},'interpreter','latex','fontsize',12,'interpreter','latex')
%
set(gca,'ticklabelinterpreter','latex','tickdir','out','fontsize',14,'ylim',[0 ceil(max(RD(hr).SUU+RD(hr).SVV))])
%
lg = legend('$S_\mathrm{uu}+S_\mathrm{vv}$','$S_\mathrm{uu}+S_\mathrm{vv}-(g/h)S_{\eta\eta}$');
set(lg,'interpreter','latex')
exportgraphics(fig4,[info.rootSim,filesep,'figures',filesep,info.rootName,'ROD_velocity_spectra_obs_mod.pdf'])
%
mooringFile = [info.rootMat,filesep,info.rootName,'_virtual_moorings.mat'];
ROD_OBS = RD(hr);
save(mooringFile,'ROD_MOD','ROD_OBS')
info.mooringFile = mooringFile;
save(info.fileName,'-struct','info')
