botFile = 'funwave_barRip0_h10t10s10d00_dep.nc';
dep = ncread(botFile,'dep');
x   = ncread(botFile,'x');
y   = ncread(botFile,'y');

momFile = 'funwave_barRip1_h10t10s10d00_MomentumTerms.nc';
Umean   = ncread(momFile,'umean');
Vmean   = ncread(momFile,'vmean');
ETAmean = ncread(momFile,'etamean');

Hmean = dep+ETAmean;
Umean = mean(Umean.*Hmean,3)./mean(Hmean,3);
Vmean = mean(Vmean.*Hmean,3)./mean(Hmean,3);
Hmean = mean(Hmean,3);

mask = Hmean>0.01;
% do helmholtz decomposition:
[psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(Umean.*Hmean,Vmean.*Hmean,y(2)-y(1),x(2)-x(1));
u_psi = u_psi./Hmean;
u_phi = u_phi./Hmean;

cm = cmocean('balance');

figure,

ax1 = subplot(4,3,[1 4 7]);
imagesc(x,y,mask.*Umean),caxis(ax1,[-0.5 0.5]), colormap(cm)
xline(50,'--k')
title('$\langle U\rangle$','interpreter','latex')
ylabel('$y$ [m]')
set(ax1,'ticklabelinterpreter','latex','xlim',[0 400])

ax2 = subplot(4,3,10);
plot(x,100*mean(mask.*Umean.*Hmean,1)./mean(Hmean,2),'-k')
xline(50,'--k')
set(ax2,'ticklabelinterpreter','latex','xlim',get(ax1,'xlim'),'ylim',[-0.1 3])
xlabel('$x$ [m]')
ylabel('avg$(u)$ [cm/s]')


ax3 = subplot(4,3,[2 5 8]);
imagesc(x,y,mask.*u_psi),caxis(ax3,[-0.5 0.5]), colormap(ax3,cm)
xline(50,'--k')
title('$\langle U\rangle_\psi$','interpreter','latex')
ylabel('$y$ [m]')
set(ax3,'ticklabelinterpreter','latex','yticklabel',[])

ax4 = subplot(4,3,11);
plot(x,100*mean(mask.*u_psi.*Hmean,1)./mean(Hmean,2),'-k')
xline(50,'--k')
set(ax4,'ticklabelinterpreter','latex','xlim',get(ax1,'xlim'),'ylim',get(ax2,'ylim'))
xlabel('$x$ [m]')

ax5 = subplot(4,3,[3 6 9]);
imagesc(x,y,mask.*u_phi),caxis(ax5,[-0.5 0.5]), colormap(ax5,cm)
xline(50,'--k')
title('$\langle U\rangle_\phi$','interpreter','latex')
ylabel('$y$ [m]')
set(ax5,'ticklabelinterpreter','latex','yticklabel',[])

ax6 = subplot(4,3,12);
plot(x,100*mean(mask.*u_phi.*Hmean,1)./mean(Hmean,2),'-k')
xline(50,'--k')
set(ax6,'ticklabelinterpreter','latex','xlim',get(ax1,'xlim'),'ylim',get(ax2,'ylim'))
xlabel('$x$ [m]')

exportgraphics(gcf,'../figures/Umean_Upsi_Uphi.pdf')