function fout = plot_FUNWAVE_run_momentum(info)
%
% USAGE: fout = plot_FUNWAVE_run_momentum(info)
%
% graphically display momentum terms:
% 1) 2D rms/mean fields of all calculated terms
% 2) 2D mean fields (velocity, vorticity)
% 3) 1D rms/mean lines of domanant cross- and along-shore terms
% 4) 1D rms/mean lines of domanant along-shore terms
%    calculated at inner- (Wsz/3), mid-SZ (Wsz/2), and at breakpoint (Wsz).
figDIR = [info.rootMOD,filesep,'figures/'];
if ~exist(figDIR,'dir')
    eval(['!mkdir -p ',figDIR])
end
fout = {};
%
% momentum file
momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
%
% $$$ %
% load from the momentum file:
x   = ncread(momFile,'x');
y   = ncread(momFile,'y');
t   = ncread(momFile,'t');
%
% load the depth for 2D plots
depFile = [info.rootMat,info.rootName,'dep.nc'];
h = ncread(depFile,'dep');
%
%
%% 0) need to check for a land-mask
if ~isfield(info,'mask')
% $$$     % full domain or sub-domain?
% $$$     if ~isfield(info,'subDomain')
% $$$         iX = find(x0>=0  & x0<=400);
% $$$         sdx = length(iX);
% $$$         sdy = length(y0);
% $$$         subDomain = [1 sdy iX(1) iX(end)];
% $$$     else
% $$$         subDomain = info.subDomain;
% $$$     end
% $$$     eta = ncread(momFile,'etamean',[subDomain([1 3]) 1] , [subDomain([2 4]) inf]);
    eta = ncread(momFile,'etamean');
    info.mask = (mean(eta,3,'omitnan')+h)>=0.1;
    % if no shorline location... use the depth-mask
    if ~isfield(info,'x_shoreline')
        tmp = repmat( 1:length(x) , length(y), 1);
        [~,iSL] = min( abs( mean(eta,3,'omitnan')+h-0.1 ), [], 2);
        info.x_shoreline = x(iSL);
    end
    save(info.fileName,'-struct','info')
end
%
% 4.1) estimate Lsz as location where cumulative integral of rms(Fbr) is 90.0% of the total
tmp     = ncread(momFile,'BrkDissX');
tmp(~info.mask)=nan;
rmsF    = std(tmp,[],[1 3],'omitnan');% this is from (3) above
rmsFtot = sum   (rmsF,'omitnan');
rmsFcum = cumsum(rmsF,'omitnan');
rmsFfrac= rmsFcum./rmsFtot;
% breakpoint at 0.9;
iBP     = find(rmsFfrac>0.9,1,'first');
xBP     = x(iBP);
%
info.x_breakpoint = xBP;
save(info.fileName,'-struct','info')
%
% get shoreline location from info-file:
xSL     = mean(info.x_shoreline);
iSL     = find(x>=xSL,1,'first');
Wsz     = xBP-xSL;
%
iINN    = find(x>=xSL+Wsz/3,1,'first');
iMID    = find(x>=xSL+Wsz*2/3,1,'first');
%
%
%% 1) make a plot of mean of BrkDissX and DxSxx 
% 1.1) time average and plot 2D fields:
vars = {'PgrdX','PgrdY','BrkDissX','BrkDissY','DxSxx','DxSxy','DySxy','DySyy','DxUUH','DxUVH','DyUVH','DyVVH','FRCX','FRCY'};
lbls = {'$gH\partial_x \bar\eta$','$gH\partial_y \bar\eta$','$-H\bar{F}_\mathrm{br,x}$','$-H\bar{F}_\mathrm{br,y}$',...
        '$\partial_x S_{xx}$','$\partial_x S_{xy}$','$\partial_y S_{xy}$','$\partial_y S_{yy}$',...
        '$\partial_x (U^2 H)$','$\partial_x (UVH)$','$\partial_y (UVH)$','$\partial_y (V^2 H)$','-$\tau_x$','-$\tau_y$'};
sgn  = {1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1 1};
%
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.2;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cpos1 = [xm+pw+ag ym       2*ag ph/2];
cpos2 = [xm+pw+ag ym+ph+ag 2*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ag+2*ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('balance');
cm2   = cmocean('amp');
% limits/grading
clim1  = [-1 1];
clim2  = [0 1];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
clrs2 = clim2(1):diff(clim2)/255:clim2(2);
% loop over variables
for jj=1:length(vars)
    var = ncread(momFile,vars{jj});% here is where you'd depth average...
    tmp = mean(var,3,'omitnan')';
    rng    = ceil(log10(range(tmp(:))/2));
    scale  = 10^rng;
    %
    clf(fig)
    a1 = axes('units','centimeters','position',ppos1);
    imagesc(y,x,sgn{jj}*tmp/scale)
    hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
    colormap(a1,cm1),caxis(a1,clim1)
    xlabel('$y$ [m]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
    a2 = axes('units','centimeters','position',ppos2);
    imagesc(y,x,std(var,[],3,'omitnan')'/scale)
    hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)    
    colormap(a2,cm2),caxis(a2,clim2)    
    ylabel('$x$ [m]','interpreter','latex')
    set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    c1 = axes('units','centimeters','position',cpos1);
    imagesc(0,clrs1,reshape(cm1,256,1,3))
    xlabel(c1,{['avg(',lbls{jj},')'];sprintf('[m/s]$^2\\times 10^{%d}$~',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
    set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',6,'tickdir','out','ydir','normal')
    c2 = axes('units','centimeters','position',cpos2);
    imagesc(0,clrs2,reshape(cm2,256,1,3))
    xlabel(c2,{['std(',lbls{jj},')'];sprintf('[m/s]$^2\\times 10^{%d}$~',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
    set(c2,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',6,'tickdir','out','ydir','normal')
    figname = [figDIR,info.runName,'_xshore_momentum_term_',vars{jj},'_time_averaged.pdf'];
    exportgraphics(fig,figname)
    fout = cat(1,fout,figname);
    clear var
end
close(fig)
%
%% 2) plot the time-averaged velocity:
Umean = ncread(momFile,'umean');
Vmean = ncread(momFile,'vmean');
ETAmean = ncread(momFile,'etamean');
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.2;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cpos1 = [xm+pw+ag ym       2*ag ph/2];
cpos2 = [xm+pw+ag ym+ph+ag 2*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ag+2*ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('balance');
clim1 = [-0.5 0.5];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,x,mean(Umean,3,'omitnan')')
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,mean(Vmean,3,'omitnan')')
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,'$\langle{U}\rangle$ [m/s]','interpreter','latex','fontsize',8,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
c2 = axes('units','centimeters','position',cpos2);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c2,'$\langle{V}\rangle$ [m/s]','interpreter','latex','fontsize',8,'horizontalalignment','left')
set(c2,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_velocity.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
ps    = [2*xm+pw+6*ag  2*ym+ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('balance');
clim1 = [-0.1 0.1];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,x,mean(ETAmean,3,'omitnan')')
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,'$\langle{\eta}\rangle$ [m]','interpreter','latex','fontsize',8,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_waterlevel.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%% plot time averaged vorticity:
[dUdy,~] = gradientDG(mean(Umean,3,'omitnan')./info.dy);
[~,dVdx] = gradientDG(mean(Vmean,3,'omitnan')./info.dx);
vort = dVdx-dUdy;
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.2;
ppos1 = [xm       ym         pw ph];
cpos1 = [xm+pw+ag ym       2*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('curl');
clim1 = [-0.1 0.1];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,x,vort')
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,'$\bar{\omega}$ [1/s]','interpreter','latex','fontsize',8)
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_vorticity.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%% 3) time and alongshore average dominant cross-shore terms
BrkDissX = ncread(momFile,'BrkDissX');
DxSxx    = ncread(momFile,'DxSxx');
PgrdX    = ncread(momFile,'PgrdX');
%
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.2;
ppos1 = [xm ym       pw ph];
ppos2 = [xm ym+ph+ag pw ph];
ps    = [2*xm+pw 2*ym+ag+2*ph];
%
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
tmp0 = PgrdX;tmp0(~info.mask)=nan; avg0 = mean(tmp0,[1 3],'omitnan'); 
tmp1 = DxSxx;tmp1(~info.mask)=nan; avg1 = mean(tmp1,[1 3],'omitnan');
tmp2 =-BrkDissX;tmp2(~info.mask)=nan; avg2 = mean(tmp2,[1 3],'omitnan');
tmp3 = PgrdX+DxSxx-BrkDissX;tmp3(~info.mask)=nan; avg3 = mean(tmp3,[1 3],'omitnan'); 
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,avg0,'k',x,avg1,'b',x,avg2,'r',x,avg3,'--c','linewidth',2);
xlabel('$x$ [m]','interpreter','latex')
ylabel('mean() [m/s]$^2$','interpreter','latex')
legend(p1(1:3),{'$gH\partial_x\, \bar{\eta}$','$\partial_x\, S_{xx}$','$-\bar{F}_{\mathrm{br},x}$'},'interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
clear avg0 avg1 avg2 avg3
rms0 = std(tmp0,[],[1 3],'omitnan');
rms1 = std(tmp1,[],[1 3],'omitnan');
rms2 = std(tmp2,[],[1 3],'omitnan');
clear tmp0 tmp1 tmp2 tmp3
a2 = axes('units','centimeters','position',ppos2);
plot(x, rms0,'k',x, rms1,'b',x, rms2,'r', 'linewidth',2)
ylabel('std() [m/s]$^2$','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
figname = [figDIR,info.runName,'_xshore_momentum_terms_time_and_alongshore_averaged.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%
%% 4) Alongshore terms at locations: inner = 1/4*Lsz, mid = 1/2*Lsz, brk-pt = Lsz
%
% 4.2.1) estimate time-averages of:
%        cross-shore terms:
%             advection, 
tmp1 = ncread(momFile,'DxUUH');
tmp2 = ncread(momFile,'DyUVH');
ADX  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan'); 
ADXstd= std(tmp1,[],3,'omitnan') + std(tmp2,[],3,'omitnan'); clear tmp1 tmp2
%             pressure grad,
PGX   = mean(PgrdX,3,'omitnan');
PGXstd= std (PgrdX,[],3,'omitnan');
%             radiation stress+BrkDissX,
tmp1 = ncread(momFile,'DySxy');
RSX  = mean(tmp1,3,'omitnan') + mean(DxSxx,3,'omitnan') - mean(BrkDissX,3,'omitnan');
RSXstd = std(tmp1,[],3,'omitnan') + std(DxSxx,[],3,'omitnan') + std(BrkDissX,[],3,'omitnan');
clear tmp1
% $$$ %             friction.
% $$$ FRCX = ncread(momFile,'FRCX');
% $$$ FRCX = mean(FRCX,3,'omitnan');
%        along-shore terms:
%             advection, 
tmp1 = ncread(momFile,'DyVVH');
tmp2 = ncread(momFile,'DxUVH');
ADY  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan'); 
ADYstd= std(tmp1,[],3,'omitnan') + std(tmp2,[],3,'omitnan'); clear tmp1 tmp2
%             pressure grad,
PgrdY = ncread(momFile,'PgrdY');
PGY   = mean(PgrdY,3,'omitnan');
PGYstd= std (PgrdY,[],3,'omitnan');
%             radiation,
tmp1 = ncread(momFile,'DySyy');
tmp2 = ncread(momFile,'DxSxy');
BrkDissY = ncread(momFile,'BrkDissY');
RSY     = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan') - mean(BrkDissY,3,'omitnan');
RSYstd  = std(tmp1,[],3,'omitnan') + std(tmp2,[],3,'omitnan') + std(BrkDissY,[],3,'omitnan');
clear tmp1 tmp2
% $$$ %             friction
% $$$ FRCY = ncread(momFile,'FRCX');
% $$$ FRCY = mean(FRCX,3,'omitnan');
% 4.2.2) loop over inner/mid/brk-pt
%
% figure parameters have same margins/gaps/etc as (3) above:
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
ylims   = info.Ly/2 + [-500 500];
NAMES   = {'Inner', 'Middle', 'Outer'};
scale = 1e-3;
iter    = 0;
for ii = [iINN iMID iBP]
    iter = iter+1;
    % 4.2.3) plot 1D transects
    clf(fig)
    a1 = axes('units','centimeters','position',ppos1);
    p1 = plot(y,PGY(:,ii)/scale,'k',y,ADY(:,ii)/scale,'b',y,RSY(:,ii)/scale,'r','linewidth',2);
    xlabel('$y$ [m]','interpreter','latex')
    ylabel('[m/s]$^2\times 10^{-3}$','interpreter','latex','fontsize',9)
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string','Alongshore','fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')
    hl = legend(p1,{'Pres.','Adv.','Wave'},'interpreter','latex','fontsize',8,'location','northeast');
    hl.AutoUpdate='off';
    hl.ItemTokenSize=[10 10];
% $$$     for kk=1:3
% $$$         p2 = icon(kk).Position;
% $$$         icon(kk).Position = [0.3 p2(2) 0];
% $$$         icon(2*(kk+1)).XData=[0.05 0.2];
% $$$     end
    set(a1,'tickdir','out','ticklabelinterpreter','latex','xlim',ylims)
    % 
    a2 = axes('units','centimeters','position',ppos2);
    p2 = plot(y,PGX(:,ii)/scale,'k',y,ADX(:,ii)/scale,'b',y,RSX(:,ii)/scale,'r','linewidth',2);
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string','Cross-shore','fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
    ylabel(sprintf('[m/s]$^2\\times 10^{%d}$',log10(scale)),'interpreter','latex','fontsize',9)
    title(a2,sprintf('%s Surfzone: $x=%3.0f$ m',NAMES{iter},x(ii)),'interpreter','latex')
    set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',ylims)
    figname = [figDIR,info.runName,'dominant_momentum_terms_',NAMES{iter},'Surfzone.pdf'];
    drawnow
    exportgraphics(fig,figname)
    fout = cat(1,fout,figname);
end
close(fig)
%    
%
% figure parameters have same margins/gaps/etc as (3) above:
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
ylims   = info.Ly/2 + [-500 500];
NAMES   = {'Inner', 'Middle', 'Outer'};
scale   = 1e-3;
iter    = 0;
for ii = [iINN iMID iBP]
    iter = iter+1;
    % 4.2.3) plot 1D transects
    clf(fig)
    scale = 1e-2;    
    a1 = axes('units','centimeters','position',ppos1);
    p1 = plot(y,PGYstd(:,ii)/scale,'k',y,ADYstd(:,ii)/scale,'b',y,RSYstd(:,ii)/scale,'r','linewidth',2);
    xlabel('$y$ [m]','interpreter','latex')
    ylabel(sprintf('[m/s]$^2\\times 10^{%d}$',log10(scale)),'interpreter','latex','fontsize',9)
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string','Alongshore: std()','fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')
    hl = legend(p1,{'Pres.','Adv.','Wave'},'interpreter','latex','fontsize',8,'location','northeast');
    hl.AutoUpdate='off';
    hl.ItemTokenSize=[10 10];
% $$$     for kk=1:3
% $$$         p2 = icon(kk).Position;
% $$$         icon(kk).Position = [0.3 p2(2) 0];
% $$$         icon(2*(kk+1)).XData=[0.05 0.2];
% $$$     end
    set(a1,'tickdir','out','ticklabelinterpreter','latex','xlim',ylims)
    % 
    a2 = axes('units','centimeters','position',ppos2);
    p2 = plot(y,PGXstd(:,ii)/scale,'k',y,ADXstd(:,ii)/scale,'b',y,RSXstd(:,ii)/scale,'r','linewidth',2);
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string','Cross-shore: std()','fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
    ylabel(sprintf('[m/s]$^2\\times 10^{%d}$',log10(scale)),'interpreter','latex','fontsize',9)
    set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'xlim',ylims)
    title(a2,sprintf('%s Surfzone: $x=%3.0f$ m',NAMES{iter},x(ii)),'interpreter','latex')
    figname = [figDIR,info.runName,'dominant_momentum_terms_std_',NAMES{iter},'Surfzone.pdf'];
    drawnow
    exportgraphics(fig,figname)
    fout = cat(1,fout,figname);
end
close(fig)
%
%% 2.2) Decompose advection terms into mean and eddy:
Hmean  = h+ETAmean;
Havg   = mean(Hmean,3,'omitnan');
%
tmp = mean(Umean.*Umean.*Hmean,3,'omitnan');
DxUUHavg = 0*Havg;
DxUUHavg(:,2:end-1) = 0.5*(tmp(:,3:end)-tmp(:,1:end-2))/info.dx;
DxUUHavg(:,[1 end]) = (tmp(:,[2 end])-tmp(:,[1 end-1]))/info.dx;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DyUVHavg = 0*Havg;
DyUVHavg(2:end-1,:) = 0.5*(tmp(3:end,:)-tmp(1:end-2,:))/info.dy;
DyUVHavg([1 end],:) = (tmp([2 end],:)-tmp([1 end-1],:))/info.dy;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DyVVHavg = 0*Havg;
DyVVHavg(2:end-1,:) = 0.5*(tmp(3:end,:)-tmp(1:end-2,:))/info.dy;
DyVVHavg([1 end],:) = (tmp([2 end],:)-tmp([1 end-1],:))/info.dy;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DxUVHavg = 0*Havg;
DxUVHavg(:,2:end-1) = 0.5*(tmp(:,3:end)-tmp(:,1:end-2))/info.dy;
DxUVHavg(:,[1 end]) = (tmp(:,[2 end])-tmp(:,[1 end-1]))/info.dy;
%
ADXavg  = (DxUUHavg+DyUVHavg);
ADXeddy = ADX-ADXavg;
ADYavg  = (DyVVHavg+DxUVHavg);
ADYeddy = ADY-ADYavg;
%
%
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.2;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cpos1 = [xm+pw+ag ym       2*ag ph/2];
cpos2 = [xm+pw+ag ym+ph+ag 2*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ag+2*ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('balance');
clim1 = [-0.5 0.5];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
scale = 1e-2;
imagesc(y,x,ADXeddy'/scale)
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string',{'Eddy Advection:'; '$\partial_x \langle \overline{u^2 H} \rangle + \partial_y \langle \overline{uv H} \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,ADXavg'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Mean Advection:'; '$\partial_x \langle u \rangle^2\langle H\rangle + \partial_y \langle u\rangle\langle v \rangle\langle H\rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s]$^2$';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_vs_eddy_x_advection.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
scale = 1e-2;
cm1   = cmocean('balance');
clim1 = [-0.5 0.5];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,x,ADYeddy'/scale)
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
           'string',{'Eddy Advection:'; '$\partial_y \langle \overline{v^2 H} \rangle + \partial_x \langle \overline{uv H} \rangle$'},...
           'fitboxtotext','on','linestyle','none','interpreter','latex',...
           'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,ADYavg'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Mean Advection:'; '$\partial_y \langle v \rangle^2\langle H\rangle + \partial_x \langle u\rangle\langle v \rangle\langle H \rangle$'},...
               'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s]$^2$';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_vs_eddy_y_advection.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%% 2.4) use a Reynolds stress analogy:
Up = Umean-mean(Umean,3,'omitnan');
Vp = Vmean-mean(Vmean,3,'omitnan');
UVavg = mean( Up.*Vp, 3,'omitnan');
%
% velocity file
uFiles = dir([info.rootMat,info.rootName,'uwavg_*.nc']);
vFiles = dir([info.rootMat,info.rootName,'vwavg_*.nc']);
u = [];
v = [];
for kk=1:length(uFiles)
    utmp = ncread([uFiles(kk).folder,filesep,uFiles(kk).name],'uwavg');
    vtmp = ncread([vFiles(kk).folder,filesep,vFiles(kk).name],'vwavg');
    u = cat(3,u,utmp);
    v = cat(3,v,vtmp);
end
clear utmp vtmp
up = u-mean(Umean,3,'omitnan');
vp = v-mean(Vmean,3,'omitnan');
uvavg = mean(up.*vp,3,'omitnan');
%
ppos3 = [xm       ym+2*(ph+ag)   pw ph];
ps    = [2*xm+pw+6*ag  2*ym+ag+3*ph];
%
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
scale = 1e-2;
cm1   = cmocean('balance');
clim1 = [-1 1];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
imagesc(y,x,uvavg'/scale)
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string',{'Reynolds Stress: $\langle u'' v'' \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,UVavg'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Reynolds Stress: $\langle \bar{u} \bar{v}\rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
a3 = axes('units','centimeters','position',ppos3);
imagesc(y,x,(mean(Umean,3,'omitnan').*mean(Vmean,3,'omitnan'))'/scale)
hold(a3,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',1)
colormap(a3,cm1),caxis(a3,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a3,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
               'string',{'Reynolds Stress: $\langle u\rangle \langle v\rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s]$^2$';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_reynolds_stress_estimate.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%% 1.2) make a video of curl(BrkDiss) and curl(DxSxx...)
BrkDissY = ncread(momFile,'BrkDissY');
DySyy    = ncread(momFile,'DySyy');
%
[dyFx,~   ] = gradientDG(BrkDissX/info.dy);
[~   ,dxFy] = gradientDG(BrkDissY/info.dx);
cFbr        = dxFy - dyFx;
[dyFx,~   ] = gradientDG(DxSxx/info.dy);
[~   ,dxFy] = gradientDG(DySyy/info.dx);
cS          = dxFy - dyFx;
%
%
ps = [2*xm+pw, 2.5*ym+ph];
fig = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
tmp0 = cFbr; tmp0(~info.mask)=nan;  rms0 = std(tmp0,[],[1 3],'omitnan');
tmp1 = cS;   tmp1(~info.mask)=nan;  rms1 = std(tmp1,[],[1 3],'omitnan');
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,rms0,'k',x,rms1,'b','linewidth',2);
xlabel('$x$ [m]','interpreter','latex')
ylabel('std() [s$^{-2}$]','interpreter','latex')
legend(p1(1:2),{'curl$(\bar{F}_\mathrm{br})$','curl$(\nabla S)$'},'interpreter','latex')
a1.YAxis.Exponent = 1;
set(a1,'tickdir','out','ticklabelinterpreter','latex')
clear rms0 rms1 
clear tmp0 tmp1
figname = [figDIR,info.runName,'_std_curl_Fbr.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
rng_cFbr  = ceil(log10(3*std(cFbr(:))/2));
rng_cS    = ceil(log10(3*std(cS(:))/2));
clims_cFbr= [-0.5 0.5]*10^rng;
clims_cS  = [-0.5 0.5]*10^rng;
%
%% 2) make a video of DxSxx and curl(DxSxx)
alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
clims = [-1 1]*1e-2;
clr_map='balance';
label1 = '$\\mathrm{curl}(\\nabla S)$ ';
label2 = {'(s$^{-2}$)'};
vidName= [figDIR,info.runName,'_curl_of_radiation_stress_gradient'];
make_1panel_video_with_bathy(vidName,x,y,t,h,cS,alims,clims_cS,clr_map,label1,label2)
fout = cat(1,fout,vidName);
%
label1 = '$\\mathrm{curl}(F_\\mathrm{br})$ ';
label2 = {'(s$^{-2}$)'};
vidName= [figDIR,info.runName,'_curl_of_breaking_force'];
make_1panel_video_with_bathy(vidName,x,y,t,h,cFbr,alims,clims_cFbr,clr_map,label1,label2)
fout = cat(1,fout,vidName);
%
return
