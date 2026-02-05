function fout = plot_FUNWAVE_offline_momentum_budget(info);
%
% USAGE: fout = plot_FUNWAVE_offline_momentum_budget(info);
%
% Use offline wave averaged velocities (U+u), 
% where \hat{u} = U + u + \tilde{u}, with \tilde{u}=wave-velocity,
% to estimate momentum budget terms:
% 
% crossshore (x): d/dx( UU + uu ) + d/dy (UV + uv)
% alongshore (y): d/dx( UV + uv ) + d/dy (VV + vv)
%
% compare these to the online calculated pressure gradients:
%
% cross-shore (x): -g d/dx \eta
% along-shore (y): -g d/dy \eta

figDIR = [info.rootMOD,filesep,'figures/'];
if ~exist(figDIR,'dir')
    eval(['!mkdir -p ',figDIR])
end
fout = {};

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
ETA   = ncread(momFile,'etamean');
U     = ncread(momFile,'umean');
V     = ncread(momFile,'vmean');
%
U = mean(U,3,'omitnan');
V = mean(V,3,'omitnan');
ETA=mean(ETA,3,'omitnan');
%
% correct for Stokes transport:
H = h+ETA;
T = mean(U.*H,1);
Ustokes = T./mean(H,1);
U       = U-Ustokes;
%
g = 9.8;
[PgrdY,PgrdX] = gradientDG(g*ETA);
PgrdY = -PgrdY./info.dy;
PgrdX = -PgrdX./info.dx;
UU = U.*U;
VV = V.*V;
UV = U.*V;
%
%
% get list of velcotiy files
files  = dir([info.rootMat,info.rootName,'uwavg*.nc']);
Nf     = length(files);
%
% preallocate:
uu = 0;
uv = 0;
vv = 0;
uuh= 0;
uvh= 0;
vvh= 0;
for ii=1:Nf
    fprintf('loading from: %s \n', files(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'uwavg_%02d.nc'],ii);
    t   = ncread(fin,'t');
    nt  = length(t);
    %
    % 3) load each variable:
    vars = {'etawavg','uwavg','vwavg'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''');'])
    end
    u   = uwavg-(U+Ustokes);
    v   = vwavg-V;
    disp('fixing eta to time-mean... bug in source code')
% $$$     dep = h+etawavg-ETA;
% $$$     dep = max(dep,0);
    dep = H;
    %
    %
    uu = uu + mean(u.*u,3,'omitnan');
    uv = uv + mean(u.*v,3,'omitnan');
    vv = vv + mean(v.*v,3,'omitnan');
    uuh= uuh+ mean(u.*u.*dep,3,'omitnan');
    uvh= uvh+ mean(u.*v.*dep,3,'omitnan');
    vvh= vvh+ mean(v.*v.*dep,3,'omitnan');
end
uu=uu/Nf;
vv=vv/Nf;
uv=uv/Nf;
uuh=uuh/Nf;
uvh=uvh/Nf;
vvh=vvh/Nf;
%
% spatially smooth all fields (5m in x, 1/2 width of ripchannel in y)
if isfield(info,'lc')
    Nflty = info.lc/info.dy;
else
    Nflty = 10/info.dy;
end
if ~mod(Nflty,2), Nflty=Nflty+1;, end
Nfltx = 10/info.dx;if ~mod(Nfltx,2), Nfltx=Nfltx+1;, end
flt  = hanning(Nflty)*hanning(Nfltx)'; flt = flt./sum(flt(:));
%
uu = conv2(uu,flt,'same');
uv = conv2(uv,flt,'same');
vv = conv2(vv,flt,'same');
%
UU = conv2(UU,flt,'same');
UV = conv2(UV,flt,'same');
VV = conv2(VV,flt,'same');
PgrdY = conv2(PgrdY,flt,'same');
%
% estimate the error in assUming constant depth:
err_uuh = uu-uuh./max(H,1e-1);
err_uvh = uv-uvh./max(H,1e-1);
err_vvh = vv-uvh./max(H,1e-1);
%
% calculate eddy momentum terms
tmp = 0*uu;
tmp(:,2:end-1) = (uu(:,[3:end]) - uu(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uu(:,[2,end]) - uu(:,[1,end-1]))./(1*info.dx);

duudx = tmp;

tmp = 0*uv;
tmp(2:end-1,:) = (uv([3:end],:) - uv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (uv([2,end],:) - uv([1,end-1],:))./(1*info.dy);

duvdy = tmp;

tmp = 0*vv;
tmp(2:end-1,:) = (vv([3:end],:) - vv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (vv([2,end],:) - vv([1,end-1],:))./(1*info.dy);

dvvdy = tmp;

tmp = 0*uv;
tmp(:,2:end-1) = (uv(:,[3:end]) - uv(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uv(:,[2,end]) - uv(:,[1,end-1]))./(1*info.dx);

duvdx = tmp;
%
% $$$ tmp = 0*uu;
% $$$ tmp(:,2:end-1) = (err_uuh(:,[3:end]) - err_uuh(:,[1:end-2]))./(2*info.dx);
% $$$ tmp(:,[1 end]) = (err_uuh(:,[2,end]) - err_uuh(:,[1,end-1]))./(1*info.dx);
% $$$ 
% $$$ duudx_err = tmp;
% $$$ 
% $$$ tmp = 0*uv;
% $$$ tmp(2:end-1,:) = (err_uvh([3:end],:) - err_uvh([1:end-2],:))./(2*info.dy);
% $$$ tmp([1 end],:) = (err_uvh([2,end],:) - err_uvh([1,end-1],:))./(1*info.dy);
% $$$ 
% $$$ duvdy_err = tmp;
% $$$ 
% $$$ tmp = 0*vv;
% $$$ tmp(2:end-1,:) = (err_vvh([3:end],:) - err_vvh([1:end-2],:))./(2*info.dy);
% $$$ tmp([1 end],:) = (err_vvh([2,end],:) - err_vvh([1,end-1],:))./(1*info.dy);
% $$$ 
% $$$ dvvdy_err = tmp;
% $$$ 
% $$$ tmp = 0*uv;
% $$$ tmp(:,2:end-1) = (err_uvh(:,[3:end]) - err_uvh(:,[1:end-2]))./(2*info.dx);
% $$$ tmp(:,[1 end]) = (err_uvh(:,[2,end]) - err_uvh(:,[1,end-1]))./(1*info.dx);
% $$$ 
% $$$ duvdx_err = tmp;
%
%
% Calculate mean momentum terms:
tmp = 0*UU;
tmp(:,2:end-1) = (UU(:,[3:end]) - UU(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UU(:,[2,end]) - UU(:,[1,end-1]))./(1*info.dx);

dUUdx = tmp;

tmp = 0*UV;
tmp(2:end-1,:) = (UV([3:end],:) - UV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (UV([2,end],:) - UV([1,end-1],:))./(1*info.dy);

dUVdy = tmp;

tmp = 0*VV;
tmp(2:end-1,:) = (VV([3:end],:) - VV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (VV([2,end],:) - VV([1,end-1],:))./(1*info.dy);

dVVdy = tmp;

tmp = 0*UV;
tmp(:,2:end-1) = (UV(:,[3:end]) - UV(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UV(:,[2,end]) - UV(:,[1,end-1]))./(1*info.dx);

dUVdx = tmp;

% $$$ %
% $$$ PGY = ncread(momFile,'PgrdY')./H;
% $$$ RSY = (ncread(momFile,'DySyy') + ncread(momFile,'DxSxy'))./H;
% $$$ FRY = ncread(momFile,'FRCY')./H;
% $$$ BRY = ncread(momFile,'BrkDissY')./H;
% $$$ PGY = conv2(mean(PGY,3,'omitnan'),flt,'same');
% $$$ RSY = conv2(mean(RSY,3,'omitnan'),flt,'same');
% $$$ FRY = conv2(mean(FRY,3,'omitnan'),flt,'same');
% $$$ BRY = conv2(mean(BRY,3,'omitnan'),flt,'same');
%
ADXeddy = duudx + duvdy;
ADXmean = dUUdx + dUVdy;
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
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string',{'Eddy Advection:'; '$\partial_x \langle \overline{u^2} \rangle + \partial_y \langle \overline{uv} \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xdir','reverse')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,ADXmean'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'xdir','reverse')
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Mean Advection:'; '$\partial_x \langle u \rangle^2 + \partial_y \langle u\rangle\langle v \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s$^2$]';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_vs_eddy_x_advection_offline_calculation.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%
%
ADYeddy = dvvdy + duvdx;
ADYmean = dVVdy + dUVdx;
%
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.2;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
ppos3 = [xm       ym+2*(ph+ag)  pw ph];
cpos1 = [xm+pw+ag ym       2*ag ph/2];
cpos2 = [xm+pw+ag ym+ph+ag 2*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+2*ag+3*ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
cm1   = cmocean('balance');
clim1 = [-0.5 0.5];
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
a1 = axes('units','centimeters','position',ppos1);
scale = 1e-3;
imagesc(y,x,ADYeddy'/scale)
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string',{'Eddy Advection:'; '$\partial_y \langle \overline{v^2} \rangle + \partial_x \langle \overline{uv} \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xdir','reverse')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,ADYmean'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'xdir','reverse')
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Mean Advection:'; '$\partial_y \langle v \rangle^2 + \partial_x \langle u\rangle\langle v \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
a3 = axes('units','centimeters','position',ppos3);
imagesc(y,x,PgrdY'/scale)
hold(a3,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a3,cm1),caxis(a3,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a3,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'xdir','reverse')
    annotation('textbox','units','centimeters','position',[ppos3(1:2)+[0 0.9].*ppos3(3:4), 0.3, 0.3],...
               'string',{'Mean Pressure Gradient:'; '$-g\partial_y \langle \eta \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s$^2$]';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_vs_eddy_y_advection_offline_calculation.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%
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
imagesc(y,x,-PgrdX'/scale)
hold(a1,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a1,cm1),caxis(a1,clim1)
xlabel('$y$ [m]','interpreter','latex')
ylabel('$x$ [m]','interpreter','latex')
    annotation('textbox','units','centimeters','position',[ppos1(1:2)+[0 0.9].*ppos1(3:4), 0.3, 0.3],...
               'string',{'Cross-shore Pressure Gradient:'; '$-g\partial_x \langle \eta \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xdir','reverse')
a2 = axes('units','centimeters','position',ppos2);
imagesc(y,x,-PgrdY'/scale)
hold(a2,'on'), contour(y,x,h',[0:1:6],'-k','linewidth',0.5)
colormap(a2,cm1),caxis(a2,clim1)    
ylabel('$x$ [m]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'xdir','reverse')
    annotation('textbox','units','centimeters','position',[ppos2(1:2)+[0 0.9].*ppos2(3:4), 0.3, 0.3],...
               'string',{'Alongshore Presure Gradient:'; '$-g\partial_y \langle \eta \rangle$'},'fitboxtotext','on','linestyle','none','interpreter','latex',...
               'fontsize',8,'backgroundcolor','none')    
c1 = axes('units','centimeters','position',cpos1);
imagesc(0,clrs1,reshape(cm1,256,1,3))
xlabel(c1,{'[m/s$^2$]';sprintf('$\\times 10^{%d}$',log10(scale))},'interpreter','latex','fontsize',6,'horizontalalignment','left')
set(c1,'ticklabelinterpreter','latex','xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
% $$$ c2 = axes('units','centimeters','position',cpos2);
% $$$ imagesc(0,clrs1,reshape(cm1,256,1,3))
% $$$ xlabel(c2,'$\bar{V}$ [m/s]','interpreter','latex','fontsize',8)
% $$$ set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
figname = [figDIR,info.runName,'_time_averaged_pressure_gradient_offline_calculation.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%
%
%
%% make a scatter plot of balance between three terms... PgrdY, duudy+duvdx, 


