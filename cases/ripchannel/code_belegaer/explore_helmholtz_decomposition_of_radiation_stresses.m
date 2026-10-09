%%
clear all
close all
addpath(genpath('~/git/funwave/'))
%% Vorticity forcing from radiation stresses, Fbr, and breaking dissipation:
%
% funwave computation of Sxx = (UUmean - WWmean + 0.5*mean(eta^2))
% funwave computation of Sxy = (UVmean)
% funwave computation of Syy = (VVmean - WWmean + 0.5*mean(eta^2))
%
% where: UUmean = mean( (U_davg-U_davg_mean)^2 * (ETA - ETAmean + Depth) )
%        VVmean = mean( (V_davg-V_davg_mean)^2 * (ETA - ETAmean + Depth) )
%        UVmean = mean( (U_davg-U_davg_mean)*(V_davg-V_davg_mean) * (ETA - ETAmean + Depth) )
%        WWmean = mean( (W_surf)^2 * (ETA - ETAmean + Depth) )
%        Wsurf  = d\eta/dt = -(P_right-P_left)/dx (Q_top - Q_bottom )/dy
%
%% 0) where are we looking? archiving?
rootDIR  = '/data2/ripchannel/'
figDIR   = '/data2/ripchannel/figures_curlWaveForce/'
outDIR   = '/data2/ripchannel/mat_data/'
if ~exist(figDIR,'dir'), eval(['!mkdir ',figDIR]), end

%% 1) need a list of run-directories
%% This section follows "~/git/funwave/cases/ripchannel/code/compile_WaveAvgVelocity_stats.m"
NAMES    = {'spreadRip-barRip0','spreadRip-barRip1','spreadRip-terRip1','highRip-barRip0-s00','highRip-barRip1-s00','highRip-terRip1-s00','highRip-terRip1-s10','highRip-barRip1-s10','highRip-barRip0-s10'};
for nn = 1:length(NAMES)
    clearvars -except rootDIR figDIR outDIR NAMES nn
    NAME     = NAMES{nn}
switch NAME
  case 'uniRip-ter2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'};
    run_dirs = {'ter2D_h10t10s02d00','ter2D_h10t10s04d00','ter2D_h10t10s10d00','ter2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'uniRip-bar2D'
    runIDs   = {'uniRip','uniRip','uniRip','uniRip'};
    run_dirs = {'bar2D_h10t10s02d00','bar2D_h10t10s04d00','bar2D_h10t10s10d00','bar2D_h10t10s20d00'};
    cblbl    = {'2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'highRip-barRip0-s00'
    runIDs   = {'highRip','spreadRip','highRip'};
    run_dirs = {'barRip0_h05t10s00d00','barRip0_h10t10s00d00','barRip0_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'};   
    run_dirs = {'barRip1_h05t10s00d00','barRip1_h10t10s00d00','barRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s00'
    runIDs   = {'highRip','spreadRip','highRip'};   
    run_dirs = {'terRip1_h05t10s00d00','terRip1_h10t10s00d00','terRip1_h15t10s00d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip0-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'barRip0_h05t10s10d00','barRip0_h10t10s10d00','barRip0_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-barRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'barRip1_h05t10s10d00','barRip1_h10t10s10d00','barRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'highRip-terRip1-s10'
    runIDs   = {'highRip','spreadRip','highRip'};    
    run_dirs = {'terRip1_h05t10s10d00','terRip1_h10t10s10d00','terRip1_h15t10s10d00'};
    cblbl    = {'0.5','1','1.5'};
    cbttl    = '$H_\mathrm{s}$ [m]';
  case 'spreadRip-barRip0'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'barRip0_h10t10s00d00','barRip0_h10t10s02d00','barRip0_h10t10s04d00','barRip0_h10t10s10d00','barRip0_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-barRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'barRip1_h10t10s00d00','barRip1_h10t10s02d00','barRip1_h10t10s04d00','barRip1_h10t10s10d00','barRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
  case 'spreadRip-terRip1'
    runIDs   = {'spreadRip','spreadRip','spreadRip','spreadRip','spreadRip'};
    run_dirs = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
    cblbl    = {'0','2','4','10','20'};
    cbttl    = '$\sigma_\theta$ [$^\circ$]';
end

N     = length(runIDs);

for ii = 1:N
    % ii=1
% $$$ runID = 'spreadRip';
% $$$ runDIR= 'barRip1_h10t10s00d00';
runID  = runIDs{ii};
runDIR = run_dirs{ii};
info   = prep_belegaer_ripchannel_info(runID,runDIR);


waves = split(info.runName,'_');
waves  = split(waves{2},{'h','t','s','d'});
%
height = (str2num(waves{2})/10);
period = (str2num(waves{3}));
spread = (str2num(waves{4}));
direction = (str2num(waves{5}));
%
depFile = [info.rootMat,'funwave_',runDIR,'_dep.nc'];
rotFile = [info.rootMat,'funwave_',runDIR,'_velocity_decomposition.nc'];
momFile = [info.rootMat,'funwave_',runDIR,'_MomentumTerms.nc'];
%
% get (x,y,dep)
x   = ncread(depFile,'x');
y   = ncread(depFile,'y');
dep = ncread(depFile,'dep');
%
% get total exchange velocity:
fileInfo = ncinfo(rotFile);
variableNames = {fileInfo.Variables.Name};
if ~ismember('Urot_mean',variableNames)
    disp('calculating helmholtz decomp on mean velocity')
    Umean = mean(ncread(momFile,'umean'),3);
    Vmean = mean(ncread(momFile,'vmean'),3);
    [~,Umean,Vmean,~,~,~]=get_vel_decomposition_reGRID(Umean,Vmean,info.dx,info.dy);
else
    Umean = ncread(rotFile,'Urot_mean');
    Vmean = ncread(rotFile,'Vrot_mean');
end
U     = ncread(rotFile,'Urot');
V     = ncread(rotFile,'Vrot');
disp('using time-averaged waterlevel... bug in source code')
ETA   = mean(ncread(momFile,'etamean'),3);
%
%
% create depth mask (min-depth-resolved=0.01m, min-depth-normalize=0.1m)
H         = dep+ETA;
mask      = H>0.01;
H(~mask)  = 0;
Hmean     = max( mean(H,3), 0.1);
%
VORT_mean = curl(x,y,Umean,Vmean);
VORT      = ncread(rotFile,'VORT');
%
%% Load radiation stress terms:
DxSxx = ncread(momFile,'DxSxx');
DySyy = ncread(momFile,'DySyy');
DxSxy = ncread(momFile,'DxSxy');
DySxy = ncread(momFile,'DxSxy');
%
% normalize by depth
DxSxx = mask.*mean(DxSxx,3)./Hmean;
DySyy = mask.*mean(DySyy,3)./Hmean;
DySxy = mask.*mean(DySxy,3)./Hmean;
DxSxy = mask.*mean(DxSxy,3)./Hmean;
%
%
% $$$ [RS_psi,RSx_rot,RSy_rot,RS_phi,RSx_phi,RSy_phi] = get_vel_decomposition_reGRID(DxSxx+DySxy,DySyy+DxSxy,info.dx,info.dy);
% $$$ %
% $$$ %
% $$$ RS_vort_force = curl(x,y,RSx_rot,RSy_rot);
flt = hamming(11)*hamming(21)'; flt = flt./sum(flt);
RS_vort_force_smooth = curl(x,y,conv2(-(DxSxx),flt,'same'),conv2(-(DySyy),flt,'same'));
%
%
% $$$ RS_vort_force_smooth = conv2(mask.*RS_vort_force,flt,'same');
%
%% Load Fbr:
FbrX = ncread(momFile,'BrkDissX');
FbrY = ncread(momFile,'BrkDissY');
FbrX = mean(FbrX,3)./Hmean;
FbrY = mean(FbrY,3)./Hmean;
%
%
Fbr_vort_force_smooth = curl(x,y,conv2(FbrX,flt,'same'),conv2(FbrY,flt,'same'));
% $$$ Fbr_vort_force_smooth = conv2(mask.*Fbr_vort_force,flt,'same');
%
g = 9.8;
%% load wave height
fileInfo = ncinfo(momFile);
variableNames = {fileInfo.Variables.Name};
if ~ismember('Hsig',variableNames)
% $$$     disp('calculating helmholtz decomp on mean velocity')
% $$$     Umean = mean(ncread(momFile,'umean'),3);
% $$$     Vmean = mean(ncread(momFile,'vmean'),3);
% $$$     [~,Umean,Vmean,~,~,~]=get_vel_decomposition_reGRID(Umean,Vmean,info.dx,info.dy);
    disp('using gamma=0.45 paramterization of wave force')
    brk_mask = conv2( sqrt( FbrX.^2 + FbrY.^2 ),flt,'same' )>0;
    gamma    = 0.45;
    gamma_sz = 0;
    gamma_br = 0;
    WaveForce= conv2( sqrt(g.*dep).*gamma.^3.*brk_mask./(4*period).*(dep>0), flt,'same');
    %
    DyWaveForce = WaveForce;
    DyWaveForce(2:end-1,:)  = (WaveForce(3:end,:)   - WaveForce(1:end-2,:)  )/(2*info.dy);
    DyWaveForce([1 end],:)  = (WaveForce([2 end],:) - WaveForce([1 end-1],:))/(  info.dy);
else
    Hs = ncread(momFile,'Hsig');
    Hs = mean(Hs,3);
    if spread==0
        fltWave = hamming(25); fltWave = fltWave'/sum(fltWave);
        Hs = conv2(Hs,fltWave,'same');
    end
    %% Funwave Dispersion Relation applied to peak frequency
    % w  = sqrt(gh) * k * sqrt( (1-a1(kh)^2)/(1-a0(kh)^2) );
    % dwdk = sqrt(gh) ( sqrt( (1-a1(kh)^2)/(1-a0(kh)^2) ) + k * ... );
    a  = -0.39;
    a1 = (a + 1/3);
    k  = wavenumber_FunwaveTVD(2*pi/period,Hmean(:));
    k  = reshape(k,size(Hmean));
    tmp = (1-a1*(k.*Hmean).^2)./(1-a*(k.*Hmean).^2);
    dwdk = sqrt(g*Hmean).*( sqrt(tmp) +...
                            0.5*k.*(tmp).^(-0.5).*( (-2*a1.*k.*Hmean.^2)./(1-a*(k.*Hmean).^2) +...
                                                    (a*k.*Hmean.^2).*(1-a1*(k.*Hmean).^2).*(1-a*(k.*Hmean).^2).^(-2) ) );
    %% full dispersion relation
    disper = sqrt(g.*k.*tanh(k.*Hmean));
    dwdk_full = 0.5*(g*tanh(k.*Hmean) + g*k.*Hmean.*sech(k.*Hmean).^2)./disper;
    dwdk_full(imag(dwdk_full)~=0)=0;
    dwdk_full = real(dwdk_full);
    %% shallow water dispersion
    cg = sqrt(9.8*Hmean.*mask);
    %% Wave energy
    E  = 9.8/16.*Hs.^2;
    %% Wave energy flux
% $$$     Ecg_shallow= E.*cg;
% $$$     Ecg_funwave= E.*dwdk;
    Ecg_full= E.*dwdk_full;
    %% Non-breaking value
    Ecg_off = sum(Ecg_full.*(Hmean>6 & Hmean<8),2)./sum((Hmean>6 & Hmean<8),2);
    %% break-point at 0.9 offshore value
    [~,idx_sz] = min( abs(Ecg_full-0.9*Ecg_off),[],2);
    idx_all    = repmat(1:size(Ecg_full,2),size(Ecg_full,1),1);
    brk_mask0  = idx_all<=idx_sz;% Ecg_full<0.9*Ecg_off;
    %% grad (E cg)
    DxEcg = 0*Ecg_full;
    DxEcg(:, 2:end-1) = (Ecg_full(:,3:end)   - Ecg_full(:,1:end-2)  )/(2*info.dx);
    DxEcg(:,[1 end])  = (Ecg_full(:,[2 end]) - Ecg_full(:,[1 end-1]))/(  info.dx);
    brk_mask1  = DxEcg>=0 & brk_mask0;
    %
    %% breaking region from Fbr:
    brk_mask2 = conv2( sqrt( FbrX.^2 + FbrY.^2 ),flt,'same' )>0.01;
    %
    %% Use mask1
    idx       = repmat(1:size(FbrX,2),size(FbrX,1),1);
    idx_brk   = max(idx.*brk_mask1,[],2);
    ind_brk   = sub2ind(size(FbrX),[1:size(FbrX,1)]',idx_brk);
    %
    gamma_sz(:,:,ii) = (Hs./Hmean);
    gamma_br(:,ii) = (Hs(ind_brk)./Hmean(ind_brk));
    Xbr(:,ii) = x(idx_brk);
    %
    %% y-derivative of dissipation
    WaveForce   = conv2(3/2*(-DxEcg./max(Hmean.*dwdk_full,0.1).*mask),flt,'same');
    DyWaveForce = WaveForce;
    DyWaveForce(2:end-1,:)  = (WaveForce(3:end,:)   - WaveForce(1:end-2,:)  )/(2*info.dy);
    DyWaveForce([1 end],:)  = (WaveForce([2 end],:) - WaveForce([1 end-1],:))/(  info.dy);
end

    figure,
    ax1 = subplot(1,3,1);
    yp  = (y-info.Ly/2)./info.lc;
    imagesc(x,yp,RS_vort_force_smooth),colormap(cmocean('balance')),caxis([-0.001 0.001])
    hold on, contour(x,yp,dep,[0 1 2 3 4],'-k')
    xlabel('$x$ [m]','interpreter','latex')
    ylabel('$(y-y_y)/L_c$ [~]','interpreter','latex')
    set(ax1,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','ylim',[-10 10])
    cb1 = colorbar('location','northoutside');
    xlabel(cb1,' curl($\nabla S$) [s$^{-2}$]~~~~~~~~~~~','interpreter','latex','fontsize',10,'horizontalalignment','center')
    set(cb1,'tickdir','out','fontsize',10,'ticklength',3*get(cb1,'ticklength'),'ticklabelinterpreter','latex')
    
    ax2 = subplot(1,3,2);
    yp  = (y-info.Ly/2)./info.lc;
    imagesc(x,yp,Fbr_vort_force_smooth),colormap(cmocean('balance')),caxis([-0.001 0.001])
    hold on, contour(x,yp,dep,[0 1 2 3 4],'-k')
    xlabel('$x$ [m]','interpreter','latex')
    set(ax2,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','ylim',[-10 10])
    cb2 = colorbar('location','northoutside');
    xlabel(cb2,' curl($F_\mathrm{br}$) [s$^{-2}$]~~~~~~~~~~~','interpreter','latex','fontsize',10,'horizontalalignment','center')  
    set(cb2,'tickdir','out','fontsize',10,'ticklength',3*get(cb2,'ticklength'),'ticklabelinterpreter','latex')
    pos = get(ax1,'position');
    ax2.Position([2 4]) = pos([2 4]);
    ax1.Position([2 4]) = pos([2 4]);

    ax3 = subplot(1,3,3);
    yp  = (y-info.Ly/2)./info.lc;
    imagesc(x,yp,-DyWaveForce),colormap(cmocean('balance')),caxis([-0.001 0.001])
    hold on, contour(x,yp,dep,[0 1 2 3 4],'-k')
    xlabel('$x$ [m]','interpreter','latex')
    set(ax3,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','ylim',[-10 10])
    cb3 = colorbar('location','northoutside');
    xlabel(cb3,' $\partial_y (c_g h)^{-1} \langle\epsilon_{br}\rangle$ [s$^{-2}$]~~~~~~~~~~~~~~~~','interpreter','latex','fontsize',10,'horizontalalignment','center')    
    set(cb3,'tickdir','out','fontsize',10,'ticklength',3*get(cb3,'ticklength'),'ticklabelinterpreter','latex')
    ax3.Position([2 4]) = pos([2 4]);    
    orient landscape

    figname = [figDIR,'curl_of_wave_force_',runID,'_',runDIR,'.pdf'];
    exportgraphics(gcf,figname)
end

end