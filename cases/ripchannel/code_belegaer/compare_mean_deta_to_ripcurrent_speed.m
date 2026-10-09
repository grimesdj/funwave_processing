clear all
close all
addpath('~/git/funwave/code/')
%
%% 1) Estimate the time/space mean sea-surface deficite in rip-channel: max( E0-E_rc )
%     1.1) time average eta
%     1.2) average over y: ambient region--> 0-1000m and 2000-3000m
%                          ripchannel    --> 1500m +/- L_rc
%% 2) Also estimate above assuming a radiation stress balance: max( S0-S_rc )
%     2.1) x-integrate: dETAdx ~ - dSxxdx / (gH)
%     2.2) average ETA_offshore + above over ambient and ripchannel regions
%
%% 3) Estimate transport in this region... this is kludgy!
%     3.1) estimate transport stream function: PSI [m^3/s]
%     3.2) find max/min and take difference... net transport: T
%     3.3) estimate cross-sectional area (A) between max/min, Urip = T/A
%
%% 4) Estimate maximum time-mean velocity in the rip-channel region.
%     --> occurs at different cross-shore locations. 
%
runBATHY = 'spreadRip'
%
prefixes = {'barRip0','barRip1','terRip1'}; 
suffixes  = {'h10t10s00d00','h10t10s02d00','h10t10s04d00','h10t10s10d00','h10t10s20d00'};
mrkrs = {'o','d','s'};
lbls  = {'Barred $L_y=100$ m','Barred $L_y=50$ m','Terraced $L_y=50$ m'}
%
Ns = length(suffixes);
Uscale_log = [];
Umax_log   = [];
Xmax_log   = [];
Urip_log   = [];
Arip_log   = [];
%
fig0 = figure;
fig00= figure;
fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
for oo = 1:length(prefixes)
runIDs   = cellstr(cell2mat(cat(2,repmat(prefixes(oo),Ns,1),repmat({'_'},Ns,1),suffixes')));
% $$$     runIDs   = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
%
N = length(runIDs);
%
%
E0 = [];
E1 = [];
S0 = [];
S1 = [];
T  = [];
D  = [];
L  = [];
A  = [];
U2 = [];
sig= {};
%
fig5 = figure;
for ii=1:N
    runID = runIDs{ii};
    info = prep_local_ripchannel_info(runBATHY,runID);
    %
    %
    momFile = dir([info.rootMat,'*',runID,'*MomentumTerms.nc']);
    momFile = [momFile(1).folder,filesep,momFile(1).name]
    if ii==1
        x = ncread(momFile,'x');
        y = ncread(momFile,'y');
        depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
        depFile = [depFile(1).folder,filesep,depFile(1).name];
        h       = ncread(depFile,'dep');
    end
    %
    % Load ETA to get E0 & E1:
    ETA = ncread(momFile,'etamean');
    H   = h+ETA;
    %
    % Load (U,V) to get T:
    U  = ncread(momFile,'umean');
    V  = ncread(momFile,'vmean');
    %
    UH = U.*H;
    %
    % enforce continuity...
    UHtot = mean(UH,1);
    Utot  = UHtot./mean(max(H,0.01),1);
    %
    % correct velocities
    UH = UH-UHtot;
    U  = U-Utot;
    UHavg = mean(UH,3,'omitnan');
    VHavg = mean(V.*H,3,'omitnan');
    %
    %
    % Load Sxx & DxSxx to get S0 and S1:
    DxSxx = ncread(momFile,'DxSxx');
    Sxx   = ncread(momFile,'Sxx')  ;
    %
    % Estimate alongshore averages:
    iY0 = (y<1e3 | y>2e3);
    iY1 = (y>info.Ly/2 -  info.lc & y<info.Ly/2 +  info.lc);
    iY2 = (y>info.Ly/2 -3*info.lc & y<info.Ly/2 +3*info.lc);
    %
    e0  = mean(ETA(iY0,:,:),[1 3],'omitnan');
    e1  = mean(ETA(iY1,:,:),[1 3],'omitnan');
    %
    if ii==1
        ETA0 = mean(ETA,3,'omitnan');
    else
        dETA = mean(ETA,3,'omitnan')-ETA0;
        figure(fig00)
        subplot(N-1,1,ii-1)
        imagesc(y-info.Ly/2,x,dETA'), caxis([-0.01 0.01]), colormap(cmocean('balance'))
        colorbar
        hold on,
        yline(50, '--r' )
        if ii==2
            title(sprintf('%s',prefixes{oo}))
        end
    end
    %
    % Estimate ETA0/ETA1 using: gH DxETA ~ -DxSxx--> ETA~ETA_offshore + \int_x0^x[DxSxx/(gH)];
    g = 9.8;
    ETA_off = mean(ETA(:,end,:)      ,3,'omitnan');
    dSxx    =-mean(cumsum(DxSxx./(g*H),2,'reverse')*info.dx, 3,'omitnan');
    %    dSxx    = mean( (Sxx - Sxx(:,end,:))./(g*H),3,'omitnan');
    s0  = mean(ETA_off(iY0,:) - dSxx(iY0,:), 1,'omitnan');% mean(DxSxx(iY0,:,:),[1 3],'omitnan');
    s1  = mean(ETA_off(iY1,:) - dSxx(iY1,:), 1,'omitnan');% mean(DxSxx(iY1,:,:),[1 3],'omitnan');
% $$$ 
% $$$     % What is the x-dependent correction to s0 and e0?
% $$$     coeff = [0*s0(x>75)'+1 s0(x>75)'] \ e0(x>75)';
% $$$     k0    = coeff(2)
% $$$     coeff = [0*s1(x>75)'+1 s1(x>75)'] \ e1(x>75)';
% $$$     k1    = coeff(2)
    %
% $$$     [PSI,UHrot,VHrot,~,~,~]=get_vel_decomposition_reGRID(UHavg,VHavg,info.dx,info.dy);
% $$$     [maxPSI,imax] = max( PSI(iY2,:), [], 1);
% $$$     [minPSI,imin] = min( PSI(iY2,:), [], 1);    
% $$$     t2  =  -sign(imax-imin).*(maxPSI-minPSI);
% $$$     %
% $$$     subplot(N,1,ii)
% $$$     imagesc(x,y(iY2),PSI(iY2,:)), caxis([-20 20]), colormap(cmocean('balance'))
    %
% $$$     tmp = y(iY2);
% $$$     l2  = tmp(max(imax,imin))-tmp(min(imax,imin));
% $$$     tmp = mean(H(iY2,:,:),3,'omitnan');
% $$$     a2  = [];
% $$$     d2  = [];
% $$$     for jj=1:length(x);
% $$$         tmp1   = tmp( min(imax(jj),imin(jj)):max(imax(jj),imin(jj)), jj );
% $$$         a2(jj) = sum(tmp1)*info.dy;
% $$$         d2(jj) = mean(tmp1);
% $$$     end
    %
    %
% $$$     % crude velocity scale
% $$$     UHrc =    UH(iY2,:,:);
% $$$     Hrc  = max(H(iY2,:,:),0.01);
% $$$     UHrc(:,x<75,:)=0;
% $$$     % crude transport scale:
% $$$     ipos = UHrc>0;
% $$$     trc  = mean(sum(UHrc.*ipos,1),3);
% $$$     arc  = mean(sum(Hrc.*ipos,1),3);
% $$$     drc  = mean(sum(Hrc.*ipos,1)./sum(ipos,1),3);
% $$$     lrc  = mean(sum(ipos,1)*info.dy,3);
% $$$     urc  = Trc./arc;
    %
    %
    % get the maximum velocity in RC:
    % first, get mean velocity profile...
    tmp0       = mean(U(iY1,:,:),[1 3]);
    ix0        = find( [tmp0(2:end)-tmp0(1:end-1)].*(x(1:end-1)'>75) > 0, 1, 'first');
% $$$     [u2,imaxU] = max(tmp0, [], [1 3],'linear');
% $$$     [maxU2,imaxU2] = max(u2);
% $$$     [r,c,l] = ind2sub(size(tmp0),imaxU(imaxU2));
% $$$     idxY2   = find(iY2);
% $$$     idxX    = find(iX);    
% $$$     yUmax      = y(idxY2(r));
% $$$     xUmax      = x(idxX(c));
    %
    % define a mask for the channel:
    mask = x'>x(ix0) & iY2;
    % get all indices of positive velocities
    u2 = 0;
    t2 = 0;
    a2 = 0;
    d2 = 0;
    l2 = 0;
    for jj=1:size(UH,3)
        tmp0 = U (:,:,jj);
        tmp1 = UH(:,:,jj);
        tmp2 = H (:,:,jj);
        ipos = tmp0>0 & mask;
        %
        % maximum transport velocity in region
        [~,imax1] = max(tmp1(:).*ipos(:));
        %
        % look for continuous region surrounding rip current.
        [bw,bl] = bwboundaries(ipos);
        % logical array for all positive velocities in connected region
        % with maximum ripcurrent speed
        in_bw   = bl==bl(imax1);
        %
        % find the maximum velcoity in region
        [~,imax0] = max(tmp0(:).*in_bw(:));
        [u2,itmp] = max([u2,tmp0(imax0)]);
        [r,c] = ind2sub(size(tmp0),imax0);
        if itmp  == 2
            iXmax   = c;
            iYmax   = r;            
        end
        %
        % average transport statistics
        t2 = max(t2,sum(tmp1(:,c).*in_bw(:,c))); %t2+sum(tmp1.*in_bw,1);
        a2 = max(a2,sum(tmp2(:,c).*in_bw(:,c)));%a2+sum(tmp2.*in_bw,1);
        d2 = max(d2,sum(tmp2(:,c).*in_bw(:,c))./sum(in_bw(:,c)));%d2+sum(tmp2.*in_bw,1)./sum(in_bw,1);
        l2 = max(l2,sum(in_bw(:,c),1)*info.dy);%l2+sum(in_bw,1)*info.dy;
        if jj==1;
           figure(fig5)
           subplot(N,1,ii)
           imagesc(x,y(iY2)-info.Ly/2,tmp0(iY2,:)), caxis([-0.5 0.5]), colormap(cmocean('balance'))
           colorbar
           hold on,
           plot( x(bw{bl(imax1)}(:,2)), y(bw{bl(imax1)}(:,1))-info.Ly/2, '--r',x(c),y(r)-info.Ly/2,'*k')
           if ii==1
               title(sprintf('%s',prefixes{oo}))
           end
        end
    end
% $$$     % average the transports, areas, lengths, etc...
% $$$     t2 = t2./size(UH,3);
% $$$     a2 = a2./size(UH,3);
% $$$     d2 = d2./size(UH,3);
% $$$     l2 = l2./size(UH,3);    
    %
    %
    E0(:,ii) = e0';
    E1(:,ii) = e1';
    S0(:,ii) = s0';
    S1(:,ii) = s1';
    T (:,ii) = t2';
    D (:,ii) = d2';
    L (:,ii) = l2';
    U2(:,ii) = u2';
    A (:,ii) = a2';
    X (:,ii) = x(iXmax);
    Y (:,ii) = y(iYmax);
    %
    % get value of spread:
    str = split(runID,'_');
    str = split(str{2},{'h','t','s','d'});
    sig{ii} = str{4};
end
figure(fig5);
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripcurrent_speed_',prefixes{oo},'.pdf'];
exportgraphics(fig5,figname)

figure(fig00);
ylabel('$x$ [m]','interpreter','latex')
xlabel('$y$ [m]','interpreter','latex')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'sealevel_difference_from_zero_spread_',prefixes{oo},'.pdf'];
exportgraphics(fig00,figname)


cm = cmocean('thermal',N+1);
cm = cm(1:N,:);

dE    = E1-E0;
ylims = round(1.1*[min(dE(x>50,:),[],'all'), max(dE(x>50,:),[],'all')]*1e3)/1e3;
fig6 = figure;
colororder(cm);
plot(x,dE,'-')
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\Delta\overline{\eta_y}$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',ylims)
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
title(sprintf('%s',prefixes{oo}))
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'sealevel_difference_from_average_eta_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)

dS    = S1-S0;
ylims = round(1.1*[min(dS(x>50,:),[],'all'), max(dS(x>50,:),[],'all')]*1e3)/1e3;

fig7 = figure;
colororder(cm);
plot(x,dS,'-')
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\Delta\tilde{\eta_y}$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','ylim',ylims)
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
title(sprintf('%s',prefixes{oo}))
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'sealevel_difference_from_radiation_stress_',prefixes{oo},'.pdf'];
exportgraphics(fig7,figname)
close(fig7)



% convert run-name to spread
spread = str2num(cell2mat(sig'));
% find maximum x-transport:
% $$$ tmp = T;
% $$$ tmp(x<75) = -999;
% $$$ [maxT,imaxT] = max(tmp,[],1);
% $$$ Xmax = x(imaxT); 
% $$$ % get corresponding array location for area/depth,
% $$$ idx  = sub2ind(size(D), imaxT,1:N);
% $$$ Arip = A(idx);
Xmax = X;
Ymax = Y;
% rip-current transport velocity scale
Arip = A;
Urip = T./A;%maxT./Arip;

% get fastest offshore speed
Umax   = U2; %max(U2,[],1);

% U-scale:
Uscale = sqrt(2*g*abs(mean(dE(x>=100 & x<=200,:))));


% $$$ ylims = round(1.1*[min(T(x>50,:),[],'all'), max(T(x>50,:),[],'all')]*1e3)/1e3;
% $$$ 
% $$$ fig8 = figure;
% $$$ colororder(cm);
% $$$ plot(x,T,'-',Xmax,T(idx),'xk')
% $$$ xline(50,'--r')
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$T_x$ [m$^3$/s]','interpreter','latex')
% $$$ set(gca,'ticklabelinterpreter','latex','ylim',ylims)
% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')
% $$$ title(sprintf('%s',prefixes{oo}))
% $$$ 
% $$$ figname = ['/data2/ripchannel/',runBATHY,filesep,'cross_shore_transport_of_ripcurrent_',prefixes{oo},'.pdf'];
% $$$ exportgraphics(fig8,figname)
% $$$ close(fig8)

% $$$ figure
% $$$ colororder(cm);
% $$$ plot(x,D,'-')
% $$$ xline(50,'--r')
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$\bar{D}$ [m]','interpreter','latex')
% $$$ set(gca,'ticklabelinterpreter','latex')
% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')
% $$$ 
% $$$ figure
% $$$ colororder(cm);
% $$$ plot(x,L,'-')
% $$$ xline(50,'--r')
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$\bar{L}$ [m]','interpreter','latex')
% $$$ set(gca,'ticklabelinterpreter','latex')
% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')

% $$$ fig9 = figure;
% $$$ colororder(cm);
% $$$ plot(x,A,'-',Xmax,Arip,'x')
% $$$ xline(50,'--r')
% $$$ xlabel('$x$ [m]','interpreter','latex')
% $$$ ylabel('$\bar{A}$ [m$^2$]','interpreter','latex')
% $$$ set(gca,'ticklabelinterpreter','latex')
% $$$ colormap(cm)
% $$$ cb = colorbar; caxis([0 N])
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
% $$$ ylabel(cb,'$\sigma_\theta$','interpreter','latex')
% $$$ title(sprintf('%s',prefixes{oo}))
% $$$ 
% $$$ figname = ['/data2/ripchannel/',runBATHY,filesep,'cross_section_area_of_ripcurrent_',prefixes{oo},'.pdf'];
% $$$ exportgraphics(fig9,figname)
% $$$ close(fig9)

% $$$ figure,
% $$$ plot(spread,maxT,'o')
% $$$ xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
% $$$ ylabel('max$(T_x)$ [m$^3$/s]','interpreter','latex')
% $$$ title(sprintf('%s',prefixes{oo}))
% $$$ 
% $$$ 
% $$$ figure,
% $$$ plot(spread,Umax,'ok',spread,5*Urip,'ob',spread,Uscale,'xr')
% $$$ xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
% $$$ ylabel('[m/s]','interpreter','latex')
% $$$ title(sprintf('%s',prefixes{oo}))
% $$$ legend({'max$(U)$','5 max$(T)/A$','$\sqrt{2g \Delta\overline{\eta_y}}$'},'interpreter','latex')


figure(fig0),
plot([0 1],[0 1],'--k')
colororder(cm)
hold on,
for ii=1:N
    plot(Uscale(ii),Umax(ii),mrkrs{oo},'color',cm(ii,:),'markerfacecolor',cm(ii,:))
end

for ii=1:N
    plot(Uscale(ii),Urip(ii),mrkrs{oo},'color',cm(ii,:))
end

figure(fig1),hold on
plot(spread,  Umax./Uscale,[mrkrs{oo},'k'])
% $$$ hold on,
% $$$ plot(spread,Urip./Uscale,[mrkrs{oo},'b'])

figure(fig2),hold on
plot(spread, Xmax, [mrkrs{oo},'k'])

figure(fig3),hold on
plot(spread, A, [mrkrs{oo},'k'])

figure(fig4),hold on
plot(spread, T, [mrkrs{oo},'k'])

Umax_log(:,oo) = Umax;
Uscale_log(:,oo)=Uscale;
Urip_log(:,oo)=Urip;
Arip_log(:,oo)=A;
end

figure(fig0)
xlabel('$\sqrt{2g\Delta\eta_y}$ [m/s]','interpreter','latex')
ylabel('$U$ [m/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'max_ripcurrent_speed_vs_velocity_scale.pdf'];
exportgraphics(fig0,figname)

figure(fig1)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$U/\sqrt{2g\Delta\eta_y}$ [~]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 1.5])
legend(lbls,'interpreter','latex','autoupdate','off')
plot(spread,0*spread+1,'--k')

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ratio_of_max_ripcurrent_speed_and_velocity_scale_vs_spread.pdf'];
exportgraphics(fig1,figname)


figure(fig2)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('argmax $U(x)$ [m]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
legend(lbls,'interpreter','latex')

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'location_of_max_ripcurrent_speed_vs_spread.pdf'];
exportgraphics(fig2,figname)

figure(fig3)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$A$ [m$^2$]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
legend(lbls,'interpreter','latex')

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'cross_section_area_of_ripcurrent_transport_vs_spread.pdf'];
exportgraphics(fig3,figname)

figure(fig4)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('$T$ [m$^3$/s]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
legend(lbls,'interpreter','latex')

figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripcurrent_transport_vs_spread.pdf'];
exportgraphics(fig4,figname)

close all
save(['/data2/ripchannel/',runBATHY,'/mat_data/mean_deta_vs_mean_rip_speed.mat'])