%% make plots of ripchannel bathymetry for OSM26
figDIR = '/data2/ripchannel/figures/';
BATHYs = {'barRip1','terRip1'};
for jj=1:length(BATHYs)
    BATHY = BATHYs{jj}
    
    % load barred:
    load(['/data2/ripchannel/spreadRip/mat_data/',BATHY,'_depth.mat']);
    info = load(['/data2/ripchannel/spreadRip/mat_data/ripchannel_run_info_',BATHY,'_h10t10s10d00.mat']);
    %
    clims = [-8 1];
    clrs  = clims(1):diff(clims)/255:clims(2);
    Nwtr = ceil(255*-clims(1)/diff(clims));
    cm1 = cmocean('deep' , Nwtr);
    cm2 = cmocean('speed',256-Nwtr+2);
    cm2 = flipud(cm2(3:end,:));
    cm3 = flipud(cmocean('thermal',-clims(1)+1));
    cm  = flipud([cm2;cm1]);
    %
    figure,
    ax1 = axes;
    surf(x,y-info.Ly/2,-h,'edgecolor','none')
    xlabel('$x$ [m]','interpreter','latex','verticalalignment','bottom')
    ylabel('$y-y_0$ [m]','interpreter','latex','verticalalignment','bottom')
    zlabel('$z_b$ [m]','interpreter','latex')                        
    colormap(ax1,cm);
    caxis(ax1,clims)
    light('position',[x(end),-info.Ly/2,10])
    view(40,15)% 30.3 17.9
    set(ax1,'ydir','normal','ticklabelinterpreter','latex','xlim',[0 600],'xtick',[0 200 400 600],'ztick',[-8:2:0],'PlotBoxAspectRatio',[1.0 2 0.69],'tickdir','out','zlim',[-9.1 0],'ylim',[-1500 1500],'fontsize',17)
    set(ax1.XAxis,'TickLabelRotation',30)
    figName = [figDIR,filesep,BATHY,'_2D.png'];
    exportgraphics(gcf,figName,'resolution',1200)
    close(gcf)
end