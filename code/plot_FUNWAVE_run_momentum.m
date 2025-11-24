function fout = plot_FUNWAVE_run_momentum(info)

figDIR = [info.rootMOD,filesep,'figures/'];
if ~exist(figDIR,'dir')
    eval(['!mkdir -p ',figDIR])
end
fout = {};
%
% momentum file
momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
%% 1) make a plot of mean of BrkDissX and DxSxx 
% 1.1) BrkDiss
x   = ncread(momFile,'x');
y   = ncread(momFile,'y');
t   = ncread(momFile,'t');
%% time average and plot 2D fields:
vars = {'PgrdX','PgrdY','BrkDissX','BrkDissY','DxSxx','DxSxy','DySxy','DySyy','DxUUH','DxUVH','DyUVH','DyVVH'};
lbls = {'$gH\partial_x \bar\eta$','$gH\partial_y \bar\eta$','$-\bar{F}_\mathrm{br,x}$','$-\bar{F}_\mathrm{br,y}$',...
        '$\partial_x S_{xx}$','$\partial_x S_{xy}$','$\partial_y S_{xy}$','$\partial_y S_{yy}$',...
        '$\partial_x (U^2 H)$','$\partial_x (UVH)$','$\partial_y (UVH)$','$\partial_y (V^2 H)$'};
sgn  = {1, 1, -1, -1, 1 1 1 1 1 1 1 1};
%
% figure parameters
xm = 2.5;
ym = 2.5;
pw = 9;
ph = 2.5;
ag = 0.1;
ppos1 = [xm       ym         pw ph];
ppos2 = [xm       ym+ph+ag   pw ph];
cpos1 = [xm+pw+ag ym       5*ag ph/2];
cpos2 = [xm+pw+ag ym+ph+ag 5*ag ph/2];
ps    = [2*xm+pw+6*ag  2*ym+ag+2*ph];
fig   = figure('units','centimeters');
fig.Position(3:4) = ps;
fig.PaperSize     = ps;
fig.PaperPosition = [0 0 ps];
% colormap
clim1  = [-0.1 0.1];
clim2  = [0 0.1];
cm1   = cmocean('balance');
cm2   = cmocean('amp');
clrs1 = clim1(1):diff(clim1)/255:clim1(2);
clrs2 = clim2(1):diff(clim2)/255:clim2(2);
%
% loop over variables
for jj=1:length(vars)
    var = ncread(momFile,vars{jj});
    clf(fig)
    a1 = axes('units','centimeters','position',ppos1);
    imagesc(y,x,sgn{jj}*mean(var,3,'omitnan')')
    colormap(a1,cm1),caxis(a1,clim1)
    xlabel('$y$ [m]','interpreter','latex')
    ylabel('$x$ [m]','interpreter','latex')
    set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
    a2 = axes('units','centimeters','position',ppos2);
    imagesc(y,x,rms(var,3,'omitnan')')
    colormap(a2,cm2),caxis(a2,clim2)    
    ylabel('$x$ [m]','interpreter','latex')
    set(a2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[])
    c1 = axes('units','centimeters','position',cpos1);
    imagesc(0,clrs1,reshape(cm1,256,1,3))
    xlabel(c1,['avg(',lbls{jj},')'],'interpreter','latex','fontsize',8)
    set(c1,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
    c2 = axes('units','centimeters','position',cpos2);
    imagesc(0,clrs2,reshape(cm2,256,1,3))
    xlabel(c2,['rms(',lbls{jj},')'],'interpreter','latex','fontsize',8)
    set(c2,'xaxislocation','top','xtick',[],'yaxislocation','right','fontsize',8,'tickdir','out','ydir','normal')
    figname = [figDIR,info.runName,'_xshore_momentum_term_',vars{jj},'_time_averaged.pdf'];
    exportgraphics(fig,figname)
    fout = cat(1,fout,figname);
    clear var
end
%
%% time and alongshore average dominant cross-shore terms
BrkDissX = ncread(momFile,'BrkDissX');
DxSxx    = ncread(momFile,'DxSxx');
PgrdX    = ncread(momFile,'PgrdX');
%
xm = 2.5;
ym = 2.5;
pw = 10;
ph = 2.5;
ag = 0.1;
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
ylabel('mean() [m/s$^2$]','interpreter','latex')
legend(p1(1:3),{'$gH\partial_x\, \bar{\eta}$','$\partial_x\, S_{xx}$','$-\bar{F}_{\mathrm{br},x}$'},'interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex')
clear avg0 avg1 avg2 avg3
rms0 = rms(tmp0,[1 3],'omitnan');
rms1 = rms(tmp1,[1 3],'omitnan');
rms2 = rms(tmp2,[1 3],'omitnan');
clear tmp0 tmp1 tmp2 tmp3
a2 = axes('units','centimeters','position',ppos2);
plot(x, rms0,'k',x, rms1,'b',x, rms2,'r', 'linewidth',2)
ylabel('rms() [m/s$^2$]','interpreter','latex')
set(a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
figname = [figDIR,info.runName,'_xshore_momentum_terms_time_and_alongshore_averaged.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
clear PgrdX
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
tmp0 = cFbr; tmp0(~info.mask)=nan;  rms0 = rms(tmp0,[1 3],'omitnan');
tmp1 = cS;   tmp1(~info.mask)=nan;  rms1 = rms(tmp1,[1 3],'omitnan');
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,rms0,'k',x,rms1,'b','linewidth',2);
xlabel('$x$ [m]','interpreter','latex')
ylabel('rms() [s$^{-2}$]','interpreter','latex')
legend(p1(1:2),{'curl$(\bar{F}_\mathrm{br})$','curl$(\nabla S)$'},'interpreter','latex')
a1.YAxis.Exponent = 1;
set(a1,'tickdir','out','ticklabelinterpreter','latex')
clear rms0 rms1 
clear tmp0 tmp1
figname = [figDIR,info.runName,'_rms_curl_Fbr.pdf'];
exportgraphics(fig,figname)
fout = cat(1,fout,figname);
close(fig)
%
%% 2) make a video of DxSxx and curl(DxSxx)
alims = [mean(info.x_shoreline) x(info.subDomain(end)) y(info.subDomain(1:2))'];
clims = [-1 1]*1e-2;
clr_map='balance';
label1 = '$\\mathrm{curl}(\\nabla S)$ ';
label2 = {'(s$^{-2}$)'};
vidName= [figDIR,info.runName,'_curl_of_radiation_stress_gradient'];
make_1panel_video(vidName,x,y,t,cS,alims,clims,clr_map,label1,label2)
fout = cat(1,fout,vidName);
%
label1 = '$\\mathrm{curl}(F_\\mathrm{br})$ ';
label2 = {'(s$^{-2}$)'};
vidName= [figDIR,info.runName,'_curl_of_breaking_force'];
make_1panel_video(vidName,x,y,t,cFbr,alims,clims,clr_map,label1,label2)
fout = cat(1,fout,vidName);
%
