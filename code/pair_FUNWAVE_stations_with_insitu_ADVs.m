function info = pair_FUNWAVE_stations_with_insitu_ADVs(info)


figDir = [info.rootSim,filesep,'figures',filesep];
figRoot= [info.runName, '_'];
if ~exist(figDir,'dir')
    eval(['!mkdir -p ', figDir])
end
stationsFile = [info.rootMat,info.rootName,'gauges.mat'];
load(stationsFile)
%
%
xs = cell2mat({stations.x});
ys = cell2mat({stations.y});
[xs0,ys0] = rotate_FRF_coords(xs,ys,-info.angle);
%
root      = '/home/derek/projects/ShortCrests/ROD/data/';
% two sets of files to load:
% 1) paros sensors
PAROSfiles     = dir([root,'AllParos',filesep,info.dateTime,'.*']);
% 2) adv sensors
ADVfiles     = dir([root,info.dateTime(1:4),filesep,info.dateTime,'.*']);
filePaths = [{PAROSfiles.folder}';{ADVfiles.folder}'];
fileNames = [{PAROSfiles.name}';{ADVfiles.name}'];
filePaths{end+1} = root;
fileNames{end+1} = ['RD_2013',info.dateTime(1:4),'.mat'];
%
% extract instrument identifiers from filename
fileNames = char(fileNames);
dum     = char(split(string(fileNames),'.'));
fileIDs = strtrim(dum(:,:,2));
%
% get instrument locations
xyfin = [root,'RODSEX_sensors_all-steve.xyz'];
fid = fopen(xyfin,'r');
dat = textscan(fid,['%s',repmat('%f',1,3)]);
instID = char(dat{1});
x      = dat{2};
y      = dat{3};
%
load('/home/derek/projects/ShortCrests/ROD/mat_data/FRF_instrument_locations_2013.mat','rodX','rodY')
%
x = [x;rodX];
y = [y;rodY];
instID=[instID;'000'];
%
% rotate to the grid angle
% find the most relevant instruments
[xa,ya] = rotate_FRF_coords(x,y,info.angle);
%
% get the instID numbers to process
[dx,ix] = min(sqrt( (xa-xs).^2 + (ya-ys).^2) ,[],1);
%
% instrument location numbers to process
stationX  = xa(ix);
stationY  = ya(ix);
stationID = instID(ix,1:2);
%
% first read DOFs file, with format:
% ID
%     xFRF     yFRF   zMSL   zBOT
fileDOF = [root,info.dateTime,'.dof'];
fid     = fopen(fileDOF);
iter    = 0;
while ~feof(fid)
    iter = iter+1;
    adv_dof_ID(iter,:)   = deblank(fgetl(fid));
    adv_dof_xyzZ (iter,:)   = str2num(fgetl(fid));
end
%
% create bare structure array to mirror stations
adv = struct([]);
% loop over files, load, save all to structure, process later
for jj = 1:size(stationID,1)
    %
    % ROD files are different
    if ismember(stationID(jj,:),'00','rows')
        fin = strtrim([root, fileNames(end,:)]);
        load(fin)
        hourROD = str2num(info.dateTime(5:6))+1;
        %
        fr = RD(hourROD).fm';
        See = RD(hourROD).SEee;
        % 1a) clark already surface corrected for depth of instrument
        % 1b) kill HF shit
        Beta    = 0.5*(1-tanh( (fr-0.5)/0.05 ));
        See     = See.*Beta;
        % 2) integrate over wind-wave band
        I  = find(fr>=1/20 & fr<=1/4);
        df  = fr(2)-fr(1);
        m0 = nansum(See(I)*df);
        m1 = nansum(fr(I).*See(I)*df);        
        H = 4*sqrt( m0 );
        Tm= m0./m1;
        % peak frequency
        [~,Ip] = max(See,[],1);
        Tp= 1./fr(Ip);
        %
        adv(jj).u_dof = 0.78;
        adv(jj).u_z   = 1.22;
        adv(jj).p_dof = 0.40;
        adv(jj).See= See;
% $$$         adv(jj).Suu= RD(hourROD).SUU;
% $$$         adv(jj).Svv= RD(hourROD).SVV;        
% $$$         adv(jj).fr = fr;
% $$$         adv(jj).Hs = H;
% $$$         adv(jj).Tm = Tm;
% $$$         adv(jj).Tp = Tp;
        adv(jj).p  = RD(hourROD).p;
        adv(jj).u  = RD(hourROD).u;% these will need to be rotated!
        adv(jj).v  = RD(hourROD).v;
        adv(jj).x  = stationX(jj);
        adv(jj).y  = stationY(jj);
        adv(jj).ID = 'ROD';
        continue
    end
    %
    % process adv/paros files
    [instFiles] = find(ismember(fileIDs(:,2:end),stationID(jj,:),'rows'));
    if instFiles==0
        adv(jj).t  = [0:1/2:3072-0.5]';
        adv(jj).ID = stationID(jj,:);
    end
    %
    for ii = 1:length(instFiles)
        fin = [strtrim(filePaths{instFiles(ii)}),filesep, strtrim(fileNames(instFiles(ii),:))];
        % get (x,y) location
        id    = fin(end-2:end);
        pq    = id(1);
        num   = id(2:end);
        %
        fid = fopen(fin,'r');
        dat = fread(fid,'float32');
        % -1) is this velocity or pressure data?
        [~,indID] = ismember(id,adv_dof_ID,'rows');
        switch pq
          case 'u'
            adv(jj).u    = dat/100;
            adv(jj).u_z  = adv_dof_xyzZ(indID,3)/100;
            adv(jj).u_dof= (adv_dof_xyzZ(indID,3)-adv_dof_xyzZ(indID,4))/100;
% $$$             [psd,fr]=mywelch(dat/100,0.5,12,0.5);            
            continue
          case 'v'
            adv(jj).v = dat/100;
            adv(jj).u_z  = adv_dof_xyzZ(indID,3)/100;
            adv(jj).u_dof= (adv_dof_xyzZ(indID,3)-adv_dof_xyzZ(indID,4))/100;
            continue
        end
        switch pq
          case 'p' 
            [~,match] = ismember([id(end-1:end),'t'],lower(instID),'rows');
            adv(jj).x = xa(match);
            adv(jj).y = ya(match);
            adv(jj).p = dat./100;
            adv(jj).p_z  = adv_dof_xyzZ(indID,3)/100;
            adv(jj).p_dof= (adv_dof_xyzZ(indID,3)-adv_dof_xyzZ(indID,4))/100;
          case 'q'
            [~,match]    = ismember([id(end-1:end),'p'],lower(instID),'rows');
            adv(jj).x    = xa(match);
            adv(jj).y    = ya(match);
            adv(jj).p_Paros    = dat./100;
            if indID~=0
                adv(jj).p_Paros_z  = adv_dof_xyzZ(indID,3)/100;
                adv(jj).p_Paros_dof= (adv_dof_xyzZ(indID,3)-adv_dof_xyzZ(indID,4))/100;
            end
        end
        adv(jj).t  = [0:1/2:3072-0.5]';
        adv(jj).ID = num;
        adv(jj).fr = fr;
    end
end
%
% now loop over advs and estimate statistics
Nadv = length(adv);
% stuff for computing spectra
Nt   = 6144;% 3072 seconds
dt   = 0.5; % 2 Hz
fs = 1./dt;
%     number of samples in average
Na = 2^floor(log2(Nt));
%     number of samples in each ensemble
Ne = 2^10;
olap = 1/2;
chnks = (Na-Ne*olap-1)/(Ne*(1-olap));
%
RODflag = 0;
x_log = [];
h_log = [];
Hs_log= [];
Tm_log= [];
Dm_log= [];
sig_log= [];
Urot_log = [];
U_log = [];
V_log = [];
for jj = 1:length(adv);
    u      = adv(jj).u;
    if isempty(u), disp(['no data for: ',adv(jj).ID]), continue, end
    v      = adv(jj).v;
    p      = adv(jj).p;
    u_dof  = adv(jj).u_dof;
    p_dof  = adv(jj).p_dof;
    h      = mean(p)+p_dof;
    if strcmp(adv(jj).ID,'ROD')
        RODflag=RODflag+1;
        % compute spectra
        % stuff for computing spectra
        Nt_ROD   = 24576;% 3072 seconds
        dt_ROD   = 0.125; % 2 Hz
        fs = 1./dt_ROD;
        %     number of samples in average
        Na_ROD = 2^floor(log2(Nt_ROD));
        %     number of samples in each ensemble
        Ne_ROD = 2^11;
        chnks_ROD = (Na_ROD-Ne*olap-1)/(Ne*(1-olap));
        [Suu,fq]=welch_method(u(1:Na_ROD),dt_ROD,chnks_ROD,0.5);
        [Svv,fq]=welch_method(v(1:Na_ROD),dt_ROD,chnks_ROD,0.5);
        [See,fq]=welch_method(p(1:Na_ROD),dt_ROD,chnks_ROD,0.5);
        [Suv,fq]=welch_cospec(u(1:Na_ROD),v(1:Na_ROD),dt_ROD,chnks_ROD,0.5);        
        [Spu,fq]=welch_cospec(p(1:Na_ROD),u(1:Na_ROD),dt_ROD,chnks_ROD,0.5);
        [Spv,fq]=welch_cospec(p(1:Na_ROD),v(1:Na_ROD),dt_ROD,chnks_ROD,0.5);
    else
        % compute spectra
        [Suu,fq]=welch_method(u(1:Na),dt,chnks,0.5);
        [Svv,fq]=welch_method(v(1:Na),dt,chnks,0.5);
        [See,fq]=welch_method(p(1:Na),dt,chnks,0.5);
        [Suv,fq]=welch_cospec(u(1:Na),v(1:Na),dt,chnks,0.5);        
        [Spu,fq]=welch_cospec(p(1:Na),u(1:Na),dt,chnks,0.5);
        [Spv,fq]=welch_cospec(p(1:Na),v(1:Na),dt,chnks,0.5);
    end
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
        % get the water-depth
        om     = 2*pi*fq;
        lom=length(om);
        g      = 9.81;
        k      = wavenumber(om,h);
        cU2U   =              cosh(k.*h)./cosh(k.*u_dof);
        cP2eta =              cosh(k.*h)./cosh(k.*p_dof);
        cU2eta = (om./(g*k)).*cU2U;
        % kill frequencies above 1Hz?
        Beta   = 0.5*(1-tanh( (fq-0.5)/0.05 ));
        % surface transform spectra
        SeUU = Suu.*(cU2eta.*Beta).^2;
        SeVV = Svv.*(cU2eta.*Beta).^2;
        SeUV = Suv.*(cU2eta.*Beta).^2;
        SePP = See.*(cP2eta.*Beta).^2;
        SePU = Spu.*cU2eta.*cP2eta.*Beta.^2;
        SePV = Spv.*cU2eta.*cP2eta.*Beta.^2;
        %
        % 
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
        %
        Ztest=SePP./(SeUU+SeVV);
        qcFlag = mean(Ztest(I));
        %
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
        Tp = m0./m1;
        %
        % estimate rotational/irrotational current magnitude
        SUU      = Suu.*(cU2U.*Beta).^2;
        SVV      = Svv.*(cU2U.*Beta).^2;
        Stot_obs = SUU + SVV;
        Sirr_obs = min(SePP*(g/h),Stot_obs);
        Srot_obs = Stot_obs - Sirr_obs;
        %
        % get the vlf/ig band
        Iobs = find(fq>0.004 & fq<0.03);
        Urot_obs = sqrt(sum(Srot_obs(Iobs)*df));
        Rrot_obs = Stot_obs/Sirr_obs;
        %
        adv(jj).h  = h;
        adv(jj).Hs = Hs;
        adv(jj).Tm = Tm;
        adv(jj).dire= [mdir1 mdir2];
        adv(jj).spread = [mspread1 mspread2];
        adv(jj).qcFlag = qcFlag;
        adv(jj).freq   = fq;
        adv(jj).See    = SePP;
        adv(jj).Suu    = SUU;
        adv(jj).Svv    = SVV;
        adv(jj).Sirr   = Sirr_obs;
        adv(jj).Srot   = Srot_obs;
        adv(jj).Urot   = Urot_obs;
        adv(jj).Rrot   = Rrot_obs;
        %
        %
        if RODflag>1
            continue
        end
        % now make some plots
        fig0 = figure;
        semilogx(stations(jj).freq,stations(jj).See,'-r',...
                 adv(jj).freq,adv(jj).See,'-k','linewidth',2)
        title_str = sprintf('ID: %s, Hs_mod: %1.1f, Hs_obs: %1.1f, Q/C flag: %1.1f',adv(jj).ID,stations(jj).Hs,adv(jj).Hs,adv(jj).qcFlag);
        legend('Mod','Obs')
        ylabel('$S_{\eta\eta}$ [m$^2$/Hz]','interpreter','latex')        
        xlabel('$f$ [Hz]','interpreter','latex')
        title(title_str,'interpreter','latex','fontsize',15)
        set(gca,'ticklabelinterpreter','latex','tickdir','out')
        figname = [figDir,figRoot,'See_MOD_vs_ADV_',adv(jj).ID,'.png'];
        exportgraphics(fig0,figname)
        %
        % now plot velocity statistics
        fig1 = figure;
        p1 = semilogx(stations(jj).freq,stations(jj).Suu+stations(jj).Svv,'--k',...
                 stations(jj).freq,stations(jj).Sirr,'--r',...
                 stations(jj).freq,stations(jj).Srot,'--b',...                 
                 adv(jj).freq,adv(jj).Suu+adv(jj).Svv,'-k',...
                 adv(jj).freq,adv(jj).Sirr,'-r',...
                 adv(jj).freq,adv(jj).Srot,'-b','linewidth',2);
        title_str = sprintf('ID: %s, Urot_mod: %1.1f, Urot_obs: %1.1f, Q/C flag: %1.1f',adv(jj).ID,stations(jj).Urot,adv(jj).Urot,adv(jj).qcFlag);
        legend([p1(1), p1(4)],'Mod','Obs')
        ylabel('$S_{uu}$ [(m/s)$^2$/Hz]','interpreter','latex')        
        xlabel('$f$ [Hz]','interpreter','latex')
        title(title_str,'interpreter','latex','fontsize',15)
        set(gca,'ticklabelinterpreter','latex','tickdir','out')
        figname = [figDir,figRoot,'Suu_MOD_vs_ADV_',adv(jj).ID,'.png'];
        exportgraphics(fig1,figname)
        %
        % keep a log of Hs
        x_log  = [x_log; adv(jj).x];
        h_log  = [h_log ; stations(jj).h , adv(jj).h];
        Hs_log = [Hs_log; stations(jj).Hs, adv(jj).Hs];
        Tm_log = [Tm_log; stations(jj).Tm, adv(jj).Tm];
        Dm_log = [Dm_log; stations(jj).dire(2), adv(jj).dire(2)];
        sig_log = [sig_log; stations(jj).spread(2), adv(jj).spread(2)];                
        Urot_log = [Urot_log; stations(jj).Urot, adv(jj).Urot];
        U_log = [U_log; mean(stations(jj).u), mean(adv(jj).u)];
        V_log = [V_log; mean(stations(jj).v), mean(adv(jj).v)];                
end
%
fig2 = figure;
ax1 = subplot(4,1,1);
p1  = plot(ax1,x_log,Hs_log(:,1),'*r',x_log,Hs_log(:,2),'*k');
ylabel(ax1,'$H_s$ [m]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out','ylim',[0 1.6])
ax2 = subplot(4,1,2);
p2  = plot(ax2,x_log,Tm_log(:,1),'*r',x_log,Tm_log(:,2),'*k');
ylabel(ax2,'$T_m$ [s]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out','ylim',[0 11])
ax3 = subplot(4,1,3);
p3  = plot(ax3,x_log,Dm_log(:,1),'*r',x_log,Dm_log(:,2),'*k');
ylabel(ax3,'$\theta_m$ [deg]','interpreter','latex')
set(ax3,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out','ylim',[-20 20])
ax4 = subplot(4,1,4);
p4  = plot(ax4,x_log,sig_log(:,1),'*r',x_log,sig_log(:,2),'*k');
ylabel(ax4,'$\sigma_\theta$ [deg]','interpreter','latex')
set(ax4,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out','ylim',[0 40])
% $$$ ax5 = subplot(5,1,5);
% $$$ p5  = plot(ax5,x_log,Urot_log(:,1),'*r',x_log,Urot_log(:,2),'*k');
% $$$ ylabel(ax5,'$U_\mathrm{rot}$ [m/s]','interpreter','latex')
% $$$ xlabel(ax5,'$x$ [m]', 'interpreter','latex')
% $$$ set(ax5,'ticklabelinterpreter','latex','tickdir','out')
orient tall
figname = [figDir,figRoot,'bulk_wave_stats_MOD_vs_ADV.png'];
exportgraphics(fig2,figname)
%
fig3 = figure;
ax1 = subplot(3,1,1);
p1  = plot(ax1,x_log,Urot_log(:,1),'*r',x_log,Urot_log(:,2),'*k');
ylabel(ax1,'$U_\mathrm{rot}$ [m/s]','interpreter','latex')
set(ax1,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out','ylim',[0 0.5])
ax2 = subplot(3,1,2);
p2  = plot(ax2,x_log,U_log(:,1),'*r',x_log,U_log(:,2),'*k');
ylabel(ax2,'$\bar{U}$ [m/s]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','xticklabel',[],'tickdir','out')
ax3 = subplot(3,1,3);
p3  = plot(ax3,x_log,V_log(:,1),'*r',x_log,V_log(:,2),'*k');
ylabel(ax3,'$\bar{V}$ [m/s]','interpreter','latex')
xlabel(ax3,'$x$ [m]', 'interpreter','latex')
set(ax3,'ticklabelinterpreter','latex','tickdir','out')
figname = [figDir,figRoot,'bulk_velocity_stats_MOD_vs_ADV.png'];
exportgraphics(fig3,figname)

save(stationsFile,'-append','adv','x_log','h_log','Urot_log','U_log','V_log','Tm_log','Dm_log','sig_log')
return
% $$$ hold on,
% $$$ Hs = [];
% $$$ for jj = 1:nADV
% $$$     figure(fig1)
% $$$     [x0,y0] = rotate_FRF_coords(adv(jj).x,adv(jj).y,-angle);
% $$$     if ~isempty(adv(jj).Hs)
% $$$         p1 = plot(x0,adv(jj).Hs,'or');
% $$$         tmp = [stations.Hs(jj),adv(jj).Hs];
% $$$     else
% $$$         tmp = [stations.Hs(jj), nan];
% $$$     end
% $$$     Hs(jj,:) = [stations.Hs(jj),adv(jj).Hs];
% $$$     if ~isempty(adv(jj).HsT)
% $$$         [x0,y0] = rotate_FRF_coords(adv(jj).xT,adv(jj).yT,-angle);
% $$$         p2 = plot(x0,adv(jj).HsT,'db');
% $$$         Hs(jj,2) = nanmean([Hs(jj,2),adv(jj).HsT]);
% $$$     end
% $$$ end
legend([p0(1) p0(2) p2(1)],'Model','Obs')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$H_\mathrm{s}$ [m]','interpreter','latex')
set(gca,'ylim',[0 1.2*max(Hs_log(:))],'ticklabelinterpreter','latex','tickdir','out')
figname = [figDir,figRoot,'Hs_MOD_vs_ADV.png'];
exportgraphics(fig1,figname)
%
% $$$ RODflag = 0;
% $$$ for jj = 1:nADV
% $$$     if (isempty(adv(jj).See) & isempty(adv(jj).SeT)) | RODflag
% $$$         continue
% $$$     elseif strcmp(adv(jj).ID,'ROD')
% $$$         RODflag=1;
% $$$     end
% $$$     figure;
% $$$     s1 = semilogx(stations(1).freq,stations(1).See(:,jj),'-k');
% $$$     hold on,
% $$$     if ~isempty(adv(jj).See) 
% $$$         s2 = semilogx(adv(jj).fr,adv(jj).See,'-r');
% $$$         legend([s1 s2],'model','paros')
% $$$     end
% $$$     if ~isempty(adv(jj).SeT)
% $$$         s3 = semilogx(adv(jj).fr,adv(jj).SeT,'--b')
% $$$         legend([s1 s2 s3],'model','paros','triton')
% $$$     end
% $$$     title(['Location: ',adv(jj).ID],'interpreter','latex')
% $$$     xlabel('$f$ [Hz]','interpreter','latex')
% $$$     ylabel('$S_{\eta\eta}$ [m$^2$/Hz]','interpreter','latex')
% $$$     set(gca,'xlim',[4e-3 1e0],'tickdir','out','ticklabelinterpreter','latex')
% $$$     figname = sprintf([figDir,figRoot,'See_MOD_vs_ADV_%s.png'],adv(jj).ID);
% $$$ %    figname = sprintf('../figures/See_avd_mod_%s.png',adv(jj).ID);                      
% $$$     exportgraphics(gcf,figname)
% $$$ end
% $$$ %
