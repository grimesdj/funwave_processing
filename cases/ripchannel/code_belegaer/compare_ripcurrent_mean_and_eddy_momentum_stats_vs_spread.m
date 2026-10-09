clear all
close all
addpath('~/git/funwave/code/')
%
%% 1) Estimate mean (UU,UV,etc) and eddy (uu,uv,etc) advection terms
%
% crossshore (x): d/dx( UU + uu ) + d/dy (UV + uv)
% alongshore (y): d/dx( UV + uv ) + d/dy (VV + vv)
%
%% 2) Estimate and compare these to pressure gradients:
%
% cross-shore (x): -g d/dx \eta
% along-shore (y): -g d/dy \eta
%
%
runBATHY = 'spreadRip'
%
prefixes = {'barRip0','barRip1','terRip1'}; 
suffixes  = {'h10t10s00d00','h10t10s02d00','h10t10s04d00','h10t10s10d00','h10t10s20d00'};
mrkrs = {'o','d','s'};
lbls  = {'Barred $L_y=100$ m','Barred $L_y=50$ m','Terraced $L_y=50$ m'}
%
Ns = length(suffixes);
%
% $$$ fig0 = figure;
fig1  = figure;% mean & eddy advection, and pressure gradient magnitudes vs spred
fig11 = figure;% bin averages of avectino vs pgrd, w/ correlations
fig111= figure;% correlation coefficient btwn mean/eddy advection and pressure gradient vs spread
fig2  = figure;
% $$$ fig3 = figure;
% $$$ fig4 = figure;
%
out = struct([]);
for oo = 1:length(prefixes)
runIDs   = cellstr(cell2mat(cat(2,repmat(prefixes(oo),Ns,1),repmat({'_'},Ns,1),suffixes')));
% $$$     runIDs   = {'terRip1_h10t10s00d00','terRip1_h10t10s02d00','terRip1_h10t10s04d00','terRip1_h10t10s10d00','terRip1_h10t10s20d00'};
%
N = length(runIDs);
%
figDIR = ['/data2/ripchannel/',runBATHY,filesep,'figures'];
if ~exist(figDIR,'dir'), eval(['!mkdir -p ', figDIR]), end
%
dETA          = [];
PGY           = [];
RSX           = [];
PGX           = [];
FRY           = [];
FRX           = [];
ADYmean       = [];
ADYeddy       = [];
ADXmean       = [];
ADXeddy       = [];
DxUU          = [];
DyUV          = [];
Dxuu          = [];
Dyuv          = [];
DyVV          = [];
DxUV          = [];
Dxuv          = [];
COH_mean      = [];
COH_eddy      = [];
Sady_mean     = [];
Sady_eddy     = [];
%
sig   = {};
%
% $$$ fig5 = figure;
for ii=1:N
    runID = runIDs{ii};
    info = prep_local_ripchannel_info(runBATHY,runID);
    %
    %
    momFile = dir([info.rootMat,'funwave_',runID,'_MomentumTerms.nc']);
    momFile = [momFile(1).folder,filesep,momFile(1).name]
    if ii==1 | ~exist('dep','var')
        x = ncread(momFile,'x');
        y = ncread(momFile,'y');
        depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
        depFile = [depFile(1).folder,filesep,depFile(1).name];
        dep     = ncread(depFile,'dep');
    end
    %
    if ~isfield(info,'Ny')
        info.Ny = length(y)+1;
    end
    % load the "means"
    % Load ETA to get E0 & E1:
    ETA = ncread(momFile,'etamean');
    H   = dep+ETA;
    ETA = mean(ETA,3,'omitnan');
    Hmean = mean(H,3,'omitnan');
    %
    % Load (U,V) to get T:
    U  = ncread(momFile,'umean');
    V  = ncread(momFile,'vmean');
    %
    % Load bottom drag:
    frx = ncread(momFile,'FRCX');
    fry = ncread(momFile,'FRCY');
    % move bottom drag to RHS of momentum equation:
    frx = -mean(frx,3,'omitnan')./Hmean;
    fry = -mean(fry,3,'omitnan')./Hmean;
    %
    % Load radiation stresses/Fbr:
    DxSxx = ncread(momFile,'DxSxx');
    DySxy = ncread(momFile,'DySyy');
    Fbr   = ncread(momFile,'BrkDissX');
    rsx   = mean(-DxSxx-DySxy+Fbr,3,'omitnan')./Hmean;
    %
    % estimate mean transport
    UH = mean(U.*H,3);
    % enforce continuity on mean transport... anti-stokes transport?
    UHtot = mean(UH,1);
    H     = mean(H,3);
    mask  = H>=0.01;
    % convert to an anti-stokes velocity
    uStokes  = UHtot./mean(max(H,0.1),1);
    %
    % correct tranport/velocity for stokes transport
    UH = UH-UHtot;
    U  = (mean(U,3)-uStokes).*mask;
    V  = (mean(V,3)).*mask;
    %
    %% calculate mean momentum terms
    g = 9.8;
    [PgrdY,PgrdX] = gradientDG(g*ETA);
    % move gradients to RHS
    PgrdY = -PgrdY./info.dy;
    PgrdX = -PgrdX./info.dx;
    UU  = U.*U;
    VV  = V.*V;
    UV  = U.*V;
    %
    % Calculate mean momentum terms:
    tmp = 0*UU;
    tmp(:,2:end-1) = (UU(:,[3:end]) - UU(:,[1:end-2]))./(2*info.dx);
    tmp(:,[1 end]) = (UU(:,[2,end]) - UU(:,[1,end-1]))./(1*info.dx);
    %
    dUUdx = tmp;
    %
    tmp = 0*UV;
    tmp(2:end-1,:) = (UV([3:end],:) - UV([1:end-2],:))./(2*info.dy);
    tmp([1 end],:) = (UV([2,end],:) - UV([1,end-1],:))./(1*info.dy);
    %
    dUVdy = tmp;
    %
    tmp = 0*VV;
    tmp(2:end-1,:) = (VV([3:end],:) - VV([1:end-2],:))./(2*info.dy);
    tmp([1 end],:) = (VV([2,end],:) - VV([1,end-1],:))./(1*info.dy);
    %
    dVVdy = tmp;
    %
    tmp = 0*UV;
    tmp(:,2:end-1) = (UV(:,[3:end]) - UV(:,[1:end-2]))./(2*info.dx);
    tmp(:,[1 end]) = (UV(:,[2,end]) - UV(:,[1,end-1]))./(1*info.dx);
    %
    dUVdx = tmp;
    %
% $$$     % load eddy statistics
% $$$     uFiles = dir([info.rootMat,'*',runID,'*uwavg*.nc']);
% $$$     vFiles = dir([info.rootMat,'*',runID,'*vwavg*.nc']);
% $$$     eFiles = dir([info.rootMat,'*',runID,'*etawavg*.nc']);
% $$$     Nf     = length(uFiles);
% $$$     uu = 0; 
% $$$     uv = 0;
% $$$     vv = 0;     
% $$$     for jj = 1:Nf
% $$$         u  = ncread( [uFiles(jj).folder,filesep,uFiles(jj).name],'uwavg');
% $$$         v  = ncread( [vFiles(jj).folder,filesep,vFiles(jj).name],'vwavg');
% $$$         % correct for mean+stokes
% $$$         u = (u - (U+uStokes)).*mask;
% $$$         v = v.*mask;
% $$$         %
% $$$         uu = uu+mean(u.*u,3,'omitnan');
% $$$         uv = uv+mean(u.*v,3,'omitnan');
% $$$         vv = uv+mean(u.*v,3,'omitnan');                
% $$$     end
% $$$     uu = uu/Nf;
% $$$     uv = uv/Nf;
% $$$     vv = vv/Nf;
    % load eddy statistics
    uFiles = dir([info.rootMat,'*',runID,'_velocity_decomposition.nc']);
    vFiles = dir([info.rootMat,'*',runID,'_velocity_decomposition.nc']);
    eFiles = dir([info.rootMat,'*',runID,'_velocity_decomposition.nc']);
    Nf     = length(uFiles);
    uu = 0; 
    uv = 0;
    vv = 0;     
    for jj = 1:Nf
        u  = ncread( [uFiles(jj).folder,filesep,uFiles(jj).name],'Urot');
        v  = ncread( [vFiles(jj).folder,filesep,vFiles(jj).name],'Vrot');
        % correct for mean+stokes
        u = u.*mask;%(u - (U+uStokes)).*mask;
        v = v.*mask;
        %
        uu = uu+mean(u.*u,3,'omitnan');
        uv = uv+mean(u.*v,3,'omitnan');
        vv = vv+mean(v.*v,3,'omitnan');                
    end
    uu = uu/Nf;
    uv = uv/Nf;
    vv = vv/Nf;
    %
    % calculate eddy momentum terms
    tmp = 0*uu;
    tmp(:,2:end-1) = (uu(:,[3:end]) - uu(:,[1:end-2]))./(2*info.dx);
    tmp(:,[1 end]) = (uu(:,[2,end]) - uu(:,[1,end-1]))./(1*info.dx);
    %
    duudx = tmp;
    %
    tmp = 0*uv;
    tmp(2:end-1,:) = (uv([3:end],:) - uv([1:end-2],:))./(2*info.dy);
    tmp([1 end],:) = (uv([2,end],:) - uv([1,end-1],:))./(1*info.dy);
    %
    duvdy = tmp;
    %
    tmp = 0*vv;
    tmp(2:end-1,:) = (vv([3:end],:) - vv([1:end-2],:))./(2*info.dy);
    tmp([1 end],:) = (vv([2,end],:) - vv([1,end-1],:))./(1*info.dy);
    %
    dvvdy = tmp;
    %
    tmp = 0*uv;
    tmp(:,2:end-1) = (uv(:,[3:end]) - uv(:,[1:end-2]))./(2*info.dx);
    tmp(:,[1 end]) = (uv(:,[2,end]) - uv(:,[1,end-1]))./(1*info.dx);
    %
    duvdx = tmp;
    %
    %
    % Estimate alongshore/cross-shore averages:
    iY0 = (y<info.Ly/2 - 5*info.lc | y>info.Ly/2 + 5*info.lc);
    iY1 = (y>info.Ly/2 -   info.lc & y<info.Ly/2 +   info.lc);
    iY2 = (y>info.Ly/2 - 3*info.lc & y<info.Ly/2 + 3*info.lc);
    iX0 = (x>137.5 & x<162.5);
    %
    e0  = mean(ETA(iY0,:,:),[1 3],'omitnan');
    e1  = mean(ETA(iY1,:,:),[1 3],'omitnan');
    %
    %
    % get value of spread:
    str = split(runID,'_');
    str = split(str{2},{'h','t','s','d'});
    sig{ii} = str{4};
    spread       (ii)   = str2num(str{4});
    dETA         (:,ii) = e1-e0;
    PGY          (:,ii) = mean(PgrdY(:,iX0),2,'omitnan');
    ADYmean      (:,ii) = mean(dVVdy(:,iX0),2,'omitnan') + mean(dUVdx(:,iX0),2,'omitnan');
    ADYeddy      (:,ii) = mean(dvvdy(:,iX0),2,'omitnan') + mean(duvdx(:,iX0),2,'omitnan');
    FRY          (:,ii) = mean(fry  (:,iX0),2,'omitnan');
    DyVV         (:,ii) = mean(dVVdy(:,iX0),2,'omitnan');
    DxUV         (:,ii) = mean(dUVdx(:,iX0),2,'omitnan');    
    Dyvv         (:,ii) = mean(dvvdy(:,iX0),2,'omitnan');
    Dxuv         (:,ii) = mean(duvdx(:,iX0),2,'omitnan');    
    PGX          (:,ii) = mean(PgrdX(iY1,:).*mask(iY1,:),1,'omitnan');
    ADXmean      (:,ii) = mean(dUUdx(iY1,:).*mask(iY1,:),1,'omitnan') + ...
                          mean(dUVdy(iY1,:).*mask(iY1,:),1,'omitnan');
    ADXeddy      (:,ii) = mean(duudx(iY1,:).*mask(iY1,:),1,'omitnan') + ...
                          mean(duvdy(iY1,:).*mask(iY1,:),1,'omitnan');        
    FRX          (:,ii) = mean(frx  (iY1,:).*mask(iY1,:),1,'omitnan');
    RSX          (:,ii) = mean(rsx  (iY1,:).*mask(iY1,:),1,'omitnan');    
    DxUU         (:,ii) = mean(dUUdx(iY1,:).*mask(iY1,:),1,'omitnan');
    DyUV         (:,ii) = mean(dUVdy(iY1,:).*mask(iY1,:),1,'omitnan');
    Dxuu         (:,ii) = mean(duudx(iY1,:).*mask(iY1,:),1,'omitnan');
    Dyuv         (:,ii) = mean(duvdy(iY1,:).*mask(iY1,:),1,'omitnan');
    %
    % include coherence:
    [coh_mean,ky,Spgy,Sadymean,Spgy_adymean] = alongshore_coherence_estimate(info,PgrdY(:,iX0),dUVdx(:,iX0));
    [coh_eddy,ky,Spgy,Sadyeddy,Spgy_adyeddy] = alongshore_coherence_estimate(info,PgrdY(:,iX0),duvdx(:,iX0));
    Crsy_mean    (:,ii) = conj(mean(Spgy_adymean,2)).*mean(Spgy_adymean,2)./(mean(Spgy,2).*mean(Sadymean,2));
    Crsy_eddy    (:,ii) = conj(mean(Spgy_adyeddy,2)).*mean(Spgy_adyeddy,2)./(mean(Spgy,2).*mean(Sadyeddy,2));
    Srsy_mean   (:,ii) = mean(Sadymean,2);
    Srsy_eddy   (:,ii) = mean(Sadyeddy,2);    
    %
    [coh_mean,ky,Spgy,Sadymean,Spgy_adymean] = alongshore_coherence_estimate(info,PgrdY(:,iX0),dVVdy(:,iX0));
    [coh_eddy,ky,Spgy,Sadyeddy,Spgy_adyeddy] = alongshore_coherence_estimate(info,PgrdY(:,iX0),dvvdy(:,iX0));
    Cady_mean    (:,ii) = conj(mean(Spgy_adymean,2)).*mean(Spgy_adymean,2)./(mean(Spgy,2).*mean(Sadymean,2));
    Cady_eddy    (:,ii) = conj(mean(Spgy_adyeddy,2)).*mean(Spgy_adyeddy,2)./(mean(Spgy,2).*mean(Sadyeddy,2));
    Sady_mean   (:,ii) = mean(Sadymean,2);
    Sady_eddy   (:,ii) = mean(Sadyeddy,2);    
    %
end

xm = 2.5;
ym = 2.5;
ag = 0.5;
ph = 2.5;
pw = 6;
ppos1 = [xm ym             pw ph];
ppos2 = [xm ym+ph+ag       pw ph];
ppos3 = [xm ym+2*(ph+ag)   pw ph];
cbpos = [xm+pw+ag ym ag 2*ph/3];
ps1   = [2*xm+pw+3*ag 2*ym+ph];
ps2   = [2*xm+pw+3*ag 2*ym+2*ph+ag];
ps3   = [2*xm+pw+3*ag 2*ym+3*ph+2*ag];
%
cm = cmocean('thermal',N+1);
cm = cm(1:N,:);
%
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps3;
set(fig6,'papersize',ps3,'paperposition',[0 0 ps3]);
colororder(cm);
%
%
flt = hamming(round(2*info.lc/info.dy)); flt = flt/sum(flt);
myfilt = @(x) conv2(x,flt,'same');
%
scale = 1e-4;
%
a3 = axes('units','centimeters','position',ppos3);
%plot((y-info.Ly/2)/info.lc,(PGY+FRY)/scale,'-','linewidth',0.1), hold on
p3 = plot((y-info.Ly/2)/info.lc,myfilt(PGY+FRY)/scale,'-')%,y-info.Ly/2,myfilt(FRY)/scale,'--');
%ylabel('[m/s$^2$]','interpreter','latex')
title(sprintf('%s',prefixes{oo}))
set(a3,'ticklabelinterpreter','latex','xlim',7*[-1 1],'xticklabel',[],'fontsize',10,'ytick',[-12:4:12])
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ (y) Pressure Gradient + Friction:'; '$-g\partial_y \langle \eta \rangle+\tau_b$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
ylims = get(a3,'ylim');
%
a2 = axes('units','centimeters','position',ppos2);
%plot(a2,(y-info.Ly/2)/info.lc,(ADYmean)/scale,'-','linewidth',0.1); hold on,
p2 = plot(a2,(y-info.Ly/2)/info.lc,myfilt(ADYmean)/scale,'-');
ylabel(sprintf('[m/s$^2$]$\\times 10^{%d}$',log10(scale)),'interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',7*[-1 1],'xticklabel',[],'fontsize',10,'ylim',ylims,'ytick',[-12:4:12])
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Mean Advection:'; '$\partial_y \langle v \rangle^2 + \partial_x \langle u \rangle\langle v\rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
                'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
%plot(a1,(y-info.Ly/2)/info.lc,(ADYeddy)/scale,'-','linewidth',0.1); hold on
p1 = plot(a1,(y-info.Ly/2)/info.lc,myfilt(ADYeddy)/scale,'-');
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
%ylabel('[m/s$^2$]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',7*[-1 1],'fontsize',10,'ylim',ylims,'ytick',[-12:4:12])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Eddy Advection:'; '$\partial_y \langle \bar{v}^2 \rangle + \partial_x \langle \bar{u}\bar{v} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2 a3],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_pressure_gradient_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
%% split the momentum terms:
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
colororder(cm);
%
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(a2,(y-info.Ly/2)/info.lc,myfilt(DyVV)/scale,'-');
ylabel(sprintf('[m/s$^2$]$\\times 10^{%d}$',log10(scale)),'interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',7*[-1 1],'xticklabel',[],'fontsize',10,'ytick',[-6:2:6])
ylims = get(a2,'ylim');
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Mean Advection Term: $\partial_y \langle v \rangle^2$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
                'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(a1,(y-info.Ly/2)/info.lc,myfilt(Dyvv)/scale,'-');
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
%ylabel('[m/s$^2$]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',7*[-1 1],'fontsize',10,'ylim',ylims,'ytick',[-6:2:6])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Eddy Advection Term: $\partial_y \langle \bar{v}^2 \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_momemtum_terms_DyVV_Dyvv_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
%
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
colororder(cm);
%
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(a2,(y-info.Ly/2)/info.lc,myfilt(DxUV)/scale,'-');
ylabel(sprintf('[m/s$^2$]$\\times 10^{%d}$',log10(scale)),'interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',7*[-1 1],'xticklabel',[],'fontsize',10,'ytick',[-2:1:2])
ylims = get(a2,'ylim');
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Mean Advection Term: $\partial_x \langle u \rangle\langle v\rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
                'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(a1,(y-info.Ly/2)/info.lc,myfilt(Dxuv)/scale,'-');
xlabel('$(y-y_0)/L_c$ [~]','interpreter','latex')
%ylabel('[m/s$^2$]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',7*[-1 1],'fontsize',10,'ylim',ylims,'ytick',[-2:1:2])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Eddy Advection Term: $\partial_x \langle \bar{u}\bar{v} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_momentum_terms_DxUV_Dxuv_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
%% coherence between Reynolds Stress and Pressure Gradient:
%
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
colororder(cm);
%
ylims = [0 1.15];
%
a2 = axes('units','centimeters','position',ppos2);
p2 = semilogx(a2,ky,Crsy_mean,'-'); hold on, plot(a2, [ky(2) ky(end)],[1 1]*(1-0.05^(1/9)),'--r')
ylabel('$\gamma^2$ [~]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',[ky(2),ky(end)],'xticklabel',[],'fontsize',10,'ylim',ylims,'ytick',[0:0.25:1])
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Coherence: $-g\partial_y\langle\eta\rangle $ \& $\partial_x \langle u \rangle\langle v\rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
p1 = semilogx(a1,ky,Crsy_eddy,'-'); hold on, plot(a1, [ky(2) ky(end)],[1 1]*(1-0.05^(1/9)),'--r')
xlabel('$k_y$ [m$^{-1}$]','interpreter','latex')
ylabel('$\gamma^2$ [~]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',[ky(2),ky(end)],'fontsize',10,'ylim',ylims,'ytick',[0:0.25:1])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Coherence: $-g\partial_y\langle\eta\rangle$ \& $\partial_x \langle \bar{u}\bar{v} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_coherence_bwtn_PgrdY_and_DxUV_Dxuv_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
%% Coherence between advection dVVdy and pressure gradient:
%
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps2;
set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
colororder(cm);
%
ylims = [0 1.15];
%
a2 = axes('units','centimeters','position',ppos2);
p2 = semilogx(a2,ky,Cady_mean,'-'); hold on, plot(a2, [ky(2) ky(end)],[1 1]*(1-0.05^(1/9)),'--r')
ylabel('$\gamma^2$ [~]','interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',[ky(2),ky(end)],'xticklabel',[],'fontsize',10,'ylim',ylims,'ytick',[0:0.25:1])
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Coherence: $-g\partial_y\langle\eta\rangle $ \& $\partial_y \langle v \rangle\langle v\rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
p1 = semilogx(a1,ky,Cady_eddy,'-'); hold on, plot(a1, [ky(2) ky(end)],[1 1]*(1-0.05^(1/9)),'--r')
xlabel('$k_y$ [m$^{-1}$]','interpreter','latex')
ylabel('$\gamma^2$ [~]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',[ky(2),ky(end)],'fontsize',10,'ylim',ylims,'ytick',[0:0.25:1])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-SZ Coherence: $-g\partial_y\langle\eta\rangle$ \& $\partial_y \langle \bar{v}\bar{v} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_coherence_bwtn_PgrdY_and_DyVV_Dyvv_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
% $$$ %% Reynolds stress spectra:
% $$$ %% coherence between Reynolds Stress and Pressure Gradient:
% $$$ %
% $$$ fig6 = figure('units','centimeters');
% $$$ fig6.Position(3:4) = ps2;
% $$$ set(fig6,'papersize',ps2,'paperposition',[0 0 ps2]);
% $$$ colororder(cm);
% $$$ %
% $$$ ylims = [0 1];
% $$$ %
% $$$ a2 = axes('units','centimeters','position',ppos2);
% $$$ p2 = loglog(a2,ky,Sady_mean,'-');
% $$$ ylabel(a2,'[(m/s$^2$)$^2$/$\Delta k_y$]~~~~~~~~~~~~~~~~~','interpreter','latex')
% $$$ set(a2,'ticklabelinterpreter','latex','ylim',[1e-11 1e-4],'xticklabel',[],'fontsize',10)
% $$$ %
% $$$ annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
% $$$            'string',{'Mid-SZ Spectra: $S(\partial_x \langle u \rangle\langle v\rangle)$'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',6,'backgroundcolor','none')    
% $$$ %
% $$$ a1 = axes('units','centimeters','position',ppos1);
% $$$ p1 = loglog(a1,ky,Sady_eddy,'-'); 
% $$$ xlabel('$k_y$ [m$^{-1}$]','interpreter','latex')
% $$$ set(a1,'ticklabelinterpreter','latex','ylim',[1e-11 1e-4],'fontsize',10)
% $$$ %
% $$$ annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
% $$$            'string',{'Mid-SZ Spectra: $S(\partial_x \langle \bar{u}\bar{v} \rangle)$'},...
% $$$            'fitboxtotext','on','linestyle','none','interpreter','latex',...
% $$$            'fontsize',6,'backgroundcolor','none')    
% $$$ %
% $$$ grid([a1 a2],'on')
% $$$ %cb = colorbar; caxis([0 N])
% $$$ cb = axes('units','centimeters','position',cbpos);
% $$$ imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
% $$$ set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
% $$$ xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
% $$$ figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_spectra_of_and_DxUV_Dxuv_',prefixes{oo},'.pdf'];
% $$$ exportgraphics(fig6,figname)
% $$$ close(fig6)
%
%% compile mean reynolds stress:
iYp = y>info.Ly/2 & y<info.Ly/2+5*info.lc;
iYm = y<info.Ly/2 & y>info.Ly/2-5*info.lc;
%
tmp_mean = myfilt(DxUV);
tmp_eddy = myfilt(Dxuv);
mag_DxUV = 0.5*(-min(tmp_mean(iYp,:)) + max(tmp_mean(iYm,:)));
mag_Dxuv = 0.5*(-min(tmp_eddy(iYp,:)) + max(tmp_eddy(iYm,:)));
%
tmp_mean = myfilt(DyVV);
tmp_eddy = myfilt(Dyvv);
mag_DyVV = 0.5*(-min(tmp_mean(iYp,:)) + max(tmp_mean(iYm,:)));
mag_Dyvv = 0.5*(-min(tmp_eddy(iYp,:)) + max(tmp_eddy(iYm,:)));
%
%
tmp_mean = myfilt(DyVV+DxUV);
tmp_eddy = myfilt(Dyvv+Dxuv);
mag_ADY_mean = 0.5*(-min(tmp_mean(iYp,:)) + max(tmp_mean(iYm,:)));
mag_ADY_eddy = 0.5*(-min(tmp_eddy(iYp,:)) + max(tmp_eddy(iYm,:)));
%
tmp_pgrd = myfilt(PGY+FRY);
mag_pgrd = 0.5*(-min(tmp_pgrd(iYp,:)) + max(tmp_pgrd(iYm,:)));
%
%
figure(fig1),hold on
plot(spread,  mag_ADY_mean,[mrkrs{oo},'b'],spread, mag_ADY_eddy,[mrkrs{oo},'r'], spread, mag_pgrd,[mrkrs{oo},'k'])
%
%
figure(fig11)
scale    = 2.5e-3;
bin_pgrd = [-1:0.1:1]*scale;
bin_ADY_eddy = [];bin2_ADY_eddy = [];
bin_ADY_mean = [];bin2_ADY_mean = [];
idx = find(iYp|iYm);
for ss = 1:size(PGY,2)
    for bb = 1:length(bin_pgrd)
        inbin = find( (tmp_pgrd(iYp|iYm,ss))>bin_pgrd(bb)-0.05*scale & (tmp_pgrd(iYp|iYm,ss))<bin_pgrd(bb)+0.05*scale );
        bin_ADY_eddy(bb,ss) = mean(tmp_eddy(idx(inbin),ss));
        bin_ADY_mean(bb,ss) = mean(tmp_mean(idx(inbin),ss));                                   
        bin2_ADY_eddy(bb,ss) = std(tmp_eddy(idx(inbin),ss));
        bin2_ADY_mean(bb,ss) = std(tmp_mean(idx(inbin),ss));                                   
    end
end
%
figure(fig11),hold on
colororder(cm), plot(bin_pgrd',bin_ADY_eddy,':s',bin_pgrd',bin_ADY_mean,'-o'),axis equal
%
%
rho_mean = diag(corr(tmp_mean(iYp|iYm,:),tmp_pgrd(iYp|iYm,:)));
rho_eddy = diag(corr(tmp_eddy(iYp|iYm,:),tmp_pgrd(iYp|iYm,:)));
figure(fig111), hold on
plot(spread,  rho_mean,[mrkrs{oo},'b'],spread, rho_eddy,[mrkrs{oo},'r'])
%
%% cross-shore momentum:
fig6 = figure('units','centimeters');
fig6.Position(3:4) = ps3;
set(fig6,'papersize',ps3,'paperposition',[0 0 ps3]);
colororder(cm);
%
ylims = [-2 2];
scale = 1e-3;
%
a3 = axes('units','centimeters','position',ppos3);
p3 = plot(x,(PGX+RSX)/scale,'-')%,x,RSX/scale,'--');
title(sprintf('%s',prefixes{oo}))
set(a3,'ticklabelinterpreter','latex','xlim',[50 400],'xticklabel',[],'fontsize',10,'ytick',[-60:10:60],'ylim',ylims)
annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
           'string',{'Mid-RC (x) Pressure Gradient + Waves:'; '$-g\partial_x \langle \eta \rangle -\partial_xS_{xx}-\partial_yS_{yy} + F_\mathrm{br,x}$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')
%
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(a2,x,ADXmean/scale,'-');
ylabel(sprintf('[m/s$^2$]$\\times 10^{%d}$',log10(scale)),'interpreter','latex')
set(a2,'ticklabelinterpreter','latex','xlim',[50 400],'xticklabel',[],'fontsize',10,'ylim',ylims,'ytick',[-6:1:6])
%
annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
           'string',{'Mid-RC Mean Advection:'; '$\partial_x \langle u \rangle^2 + \partial_y \langle u \rangle\langle v\rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(a1,x,ADXeddy/scale,'-');
xlabel('$(y-y_0)$ [m]','interpreter','latex')
%ylabel('[m/s$^2$]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',[50 400],'fontsize',10,'ylim',ylims,'ytick',[-6:1:6])
%
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Mid-RC Eddy Advection:'; '$\partial_x \langle \bar{u}^2 \rangle + \partial_y \langle \bar{u}\bar{v} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',6,'backgroundcolor','none')    
%
grid([a1 a2 a3],'on')
%cb = colorbar; caxis([0 N])
cb = axes('units','centimeters','position',cbpos);
imagesc(0,(1:N)-0.5,reshape(cm,N,1,3))
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig,'yaxislocation','right','xaxislocation','top','ydir','normal','xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',8)
xlabel(cb,'$\sigma_\theta$ [$^\circ$]','interpreter','latex','horizontalalignment','left')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_crossshore_pressure_gradient_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
%% plot sea-surface difference:
g = 9.8;
% Uscale = sqrt(2*g*abs(mean(dETA(x>=100 & x<=200,:))));
%
fig6 = figure;
colororder(cm);
p1 = plot(x,dETA,'-');
xline(50,'--r')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$\Delta \langle\eta\rangle_y$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','xlim',[50 400])
colormap(cm)
cb = colorbar; caxis([0 N])
set(cb,'ylim',[0 N],'ytick',(1:N)-0.5,'yticklabel',sig)
ylabel(cb,'$\sigma_\theta$','interpreter','latex')
title(sprintf('%s',prefixes{oo}))
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_sealevel_anomaly_dETA_',prefixes{oo},'.pdf'];
exportgraphics(fig6,figname)
close(fig6)
%
    out(oo).Name         = prefixes{oo};
    out(oo).spread       = spread;    
    out(oo).dETA         = dETA;
    out(oo).PGY          = PGY;
    out(oo).RSX          = RSX;
    out(oo).PGX          = PGX;
    out(oo).FRY          = FRY;
    out(oo).FRX          = FRX;
    out(oo).ADYmean      = ADYmean;
    out(oo).ADYeddy      = ADYeddy;
    out(oo).ADXmean      = ADXmean;
    out(oo).ADXeddy      = ADYeddy;
    out(oo).DxUU         = DxUU;
    out(oo).DyUV         = DyUV;
    out(oo).Dxuu         = Dxuu;
    out(oo).Dyuv         = Dyuv;
    out(oo).DyVV         = DyVV;
    out(oo).DxUV         = DxUV;
    out(oo).Dxuv         = Dxuv;
    out(oo).Crsy_mean     = Crsy_mean;
    out(oo).Crsy_eddy     = Crsy_eddy;
    out(oo).Srsy_mean     = Srsy_mean;
    out(oo).Srsy_eddy     = Srsy_eddy;
    out(oo).Cady_mean     = Cady_mean;
    out(oo).Cady_eddy     = Cady_eddy;
    out(oo).Sady_mean     = Sady_mean;
    out(oo).Sady_eddy     = Sady_eddy;
    out(oo).ky           = ky;
end

figure(fig1)
xlabel('$\sigma_\theta$ [$^\circ$]','interpreter','latex')
ylabel('[m/s$^2$]','interpreter','latex')
title('Maximum Reynolds Stress: $\partial_x \langle u\rangle\langle v\rangle$ (black), $\partial_x \langle\bar{u}\bar{v}\rangle (red)$','fontsize',10,'interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 1.5]*1e-4)
legend(gca().Children((2*length(prefixes)-1):-2:1),lbls,'interpreter','latex','autoupdate','off')
figname = ['/data2/ripchannel/',runBATHY,filesep,'figures',filesep,'ripchannel_alongshore_momentum_terms_DxUV_Dxuv_maximum_all.pdf'];
exportgraphics(fig1,figname)

% close all
save(['/data2/ripchannel/',runBATHY,'/mat_data/mean_vs_eddy_rip_momentum.mat'],'out')