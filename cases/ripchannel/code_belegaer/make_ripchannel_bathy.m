function info = make_ripchannel_bathy(info,plotter);
%
% usage: info = make_IDEAL_planar_bathy(info,plotter);
%
% INPUT: "info" is a structure array with fields
% rootInp:   /path/to/simulation/directory/
% rootMat:   /path/to/matlab/directory to store a .mat copy of bathy
% runName:   prefix_for_run
% dx:        grid cross-shore resolution
% dy:        grid along-shore resolution
% s:         bottom slope
% Ly:        desired alongshore domain length
% Hs:        significant wave-height
% Tp:        characteristic period
%
% additional parameters for making a Gaussian bar:
% xc:        x-location of bar from shoreline
% wc:        cross-shore width of bar
% ac:        height of bar above bottom
%
% ripchannel:
% rc:        ratio of bar height to channel depth
% lc:        width of ripchannel
%
% OUTPUT: "info" with new fields
% x:       vector of cross-shore cell centers [1 ,Nx]
% y:       vector of along-shore cell centers [Ny, 1]
% h:       depth(x,y) z_bottom=-h             [Ny,Nx]
% xWM:     cross-shore location of wave-maker
% hWM:     depth at wave-maker
% bathyFile: path to .mat of full bathy
if ~exist(info.rootSim,'dir'), eval(['!mkdir -p ', info.rootSim]), end
BATHY   = regexp(info.runName,'^.+?(?=_)','match');
BATHY   = BATHY{1};
ASCII   = [info.rootInp,filesep,BATHY,'_depth.txt'];
GAUGE   = [info.rootInp,filesep,BATHY,'_gauge.txt'];
MAT     = [info.rootMat,filesep,BATHY,'_depth.mat'];
% $$$ ASCII1D = [info.rootMOD,filesep,info.runName,'1D_depth.txt'];
% $$$ MAT1D   = [info.rootMat,filesep,info.runName,'1D_depth.mat'];
%
% get relevant vars from input structure
dx = info.dx;
dy = info.dy;
Ly = info.Ly;
s  = info.slope;
Hs = info.Hs;
Tp = info.Tp;
%
% Here we're fixing the h(x0)=0 location:
x0     = 50;
%
% depth at wave maker should satisfy k*h<=pi for all runs/frequencies, use maximum bounds
Tp_max = 10; Tp_min = 8;
Hs_max = 1.5;
hWM = 9;
% satisfy this condition
%
% for building wave-maker and sponges, use the largest period...
kWM  = wavenumber(2*pi/Tp_max,hWM);
L    = 2*pi/kWM;
delta= 2;
W    = max(delta*L/2,85);
%
% get length of sub-regions
slopeLength = ceil(hWM/s);
flatLength  = ceil(2.75*L+W);% 1L onshore of W, with ~1L offshore + ~0.75L sponge
Lx          = x0 +slopeLength + flatLength;
Nx          = round(Lx/dx)+1;
xWM         = round((Lx-1.5*L-W/2));
%
x = [0:Nx-1]*dx;
h0= (x-x0)*s;
h0(h0>=hWM)=hWM;
% $$$ % spydell & feddo 2009 bathy
% $$$ h0(h0<=0.2)=0.2;
%
% now smooth the slope near the kinks
Nw         = round(L/dx/2); if iseven(Nw),Nw=Nw+1; end
hammWindow = hamming(Nw); hammWindow = hammWindow/sum(hammWindow);
smoothKink = conv(h0,hammWindow,'same');
h0(x0+Nw:end-Nw) = smoothKink(x0+Nw:end-Nw);
%
%
% 1) add terrace/bar/channel feature
if isfield(info,'xc')
    bar = info.ac.*exp( -0.5*( (x-info.xc)/info.wc ).^2 );
    %
    % 2) if this is a 2D run, is there a channel?..
    if isfield(info,'lc')
        dy = info.dy;
        Ny = info.Ly/dy;
        y  = [0:Ny-1]'*dy;
        bar = bar - bar.*info.rc.*exp( -0.5*( (y-info.Ly/2)/info.lc ).^2);
    end
    %
    % 3) add this feature to the planar profile
    h0  = h0-bar;
end
%
% $$$ % spydell & feddo 2009 bathy
% $$$ Nw         = Nw*(2/hWM); if iseven(Nw),Nw=Nw+1; end
% $$$ hammWindow = hamming(Nw); hammWindow = hammWindow/sum(hammWindow);
% $$$ smoothKink = conv(h0,hammWindow,'same');
% $$$ h0(Nw:x0+Nw) = smoothKink(Nw:x0+Nw);
%
%
%
if exist('plotter','var')
    if plotter
        figure,
        if size(h0,1)>1
            plot(x,h0(round(Ny/2),:),'--',x,h0(1,:),'-','linewidth',2),
        else
            plot(x,h0),
        end
        set(gca,'ydir','reverse')
        ylabel('depth','interpreter','latex')
        xlabel('cross-shore','interpreter','latex')
        ylims = [min(h0(:))-0.25 hWM+0.1];
        hold on,patch( [xWM-W/2; xWM-W/2; xWM+W/2; xWM+W/2], [ylims(1); hWM; hWM; ylims(1)],'b','facealpha',0.25)
        hold on,patch( [x(end)-0.75*L; x(end)-0.75*L; x(end); x(end)], [ylims(1); hWM; hWM; ylims(1)],'g','facealpha',0.25)
        hold on,patch( [x(1)+10; x(1)+10; x(1); x(1)], [ylims(1); h0(round(10/dx)); h0(1); ylims(1)],'g','facealpha',0.25)
        set(gca,'ticklabelinterpreter','latex','fontsize',15,'xlim',[x(1) x(end)],'ylim',ylims,'plotboxaspectratio',[1 0.5 1],'tickdir','out')
        grid on
        %
        %
        figName = [info.rootInp,filesep,BATHY,'.png'];
        exportgraphics(gcf,figName)
        %
        if isfield(info,'lc')
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
            surf(x,y-info.Ly/2,-h0,'edgecolor','none')
            xlabel('$x$ [m]','interpreter','latex','verticalalignment','bottom')
            ylabel('$y-y_0$ [m]','interpreter','latex','verticalalignment','bottom')
            zlabel('$z_b$ [m]','interpreter','latex')                        
            colormap(ax1,cm);
            caxis(ax1,clims)
            light('position',[x(end),-info.Ly/2,10])
% $$$             ylabel(cb1,'$z_b$ [m]','interpreter','latex')
% $$$             set(cb1,'ticklabelinterpreter','latex')
            view(40,15)% 30.3 17.9
            set(ax1,'ydir','normal','ticklabelinterpreter','latex','xlim',[0 600],'xtick',[0 200 400 600],'ztick',[-8:2:0],'PlotBoxAspectRatio',[1.0 2 0.69],'tickdir','out','zlim',[-9.1 0],'ylim',[-600 600])
            set(ax1.XAxis,'TickLabelRotation',30)
% $$$             cb1 = colorbar;
% $$$             ax1pos = ax1.Position;
% $$$             cb1.Position(4) = 0.4*cb1.Position(4);
% $$$             ax1.Position = ax1pos;
            %
% $$$             ax2 = copyobj(ax1,gcf);
% $$$             delete(ax2.Children)
% $$$             contour(ax2,y-info.Ly/2,x,-h0',[-1:-1:-8],'-'), hold(ax2,'on')
% $$$             contour(ax2,y-info.Ly/2,x,-h0',[-0.5:-1:-8],':')
% $$$             yline(ax2,x0,'--r'),
% $$$             colormap(ax2,cm3)
% $$$             caxis(ax2,[clims(1) 0])
% $$$             set(ax2,'color','none','xcolor','none','ycolor','none','ydir','normal','position',ax1.Position,'ylim',[0 350],'plotboxaspectratio',ax1.PlotBoxAspectRatio)
% $$$             cb2 = axes;
% $$$             cb2.Position = cb1.Position;
% $$$             cb2.Position(2)  = sum([cb1.Position([2 4]), 0.1]);
% $$$             imagesc(cb2,0,[clims(1):0],reshape(cm3,size(cm3,1),1,3))
% $$$             ylabel(cb2,'$z_b$ [m]','interpreter','latex')            
% $$$             set(cb2,'ytick',[-8:2:0],'ticklabelinterpreter','latex','tickdir','out','ticklength',4*get(cb2,'ticklength'),'ylim',[clims(1)-0.5 0+0.5],'ydir','normal','yaxislocation','right','box','off','xtick',[])
% $$$             set(cb1,'ticklabelinterpreter','latex','tickdir','out','ticklength',4*get(cb1,'ticklength'))
            figName = [info.rootInp,filesep,BATHY,'_2D','.png']
            exportgraphics(gcf,figName)
        end
    end
end
%
%
info.bathyFile = MAT;
fprintf('Updating the following fields in file created by "run_info_IDEAL.m" \n')
fprintf(' info.bathyFile= %s \n', MAT)
if isfield(info,'is1D')
    h      = [h0;h0;h0];
    y      = [-dy;0;dy];
    Ny     = length(y);
    info.Ly= range(y);
    fprintf(' info.Ly        = %f \n', info.Ly)    
elseif size(h0,1)==1
    Ny = round(Ly/dy);
    y = [0:Ny-1]*dy;% the Ny^th point is identical to the 1st point
    h = repmat(h0,Ny,1);
else
    h = h0;
end
%
save(MAT,'x','y','h')
eval(['save -ASCII ',ASCII,' h'])
%
info.x = x;
info.y = y;
info.Nx = Nx;
info.Ny = Ny;
info.h = h;
info.xWM=xWM;
info.hWM=hWM;
%
fprintf(' info.Nx        = %u \n', info.Nx)
fprintf(' info.Ny        = %u \n', info.Ny)
fprintf(' info.xWM       = %f \n', info.xWM)
fprintf(' info.xWM       = %f \n', info.hWM)
%
% create gauge locations... 1D vs 2D rip-channel
ix     = find(h0(1,:)>0.2,1,'first');
xg     =[ix + 0:round(20/dx):round((x0+slopeLength)/dx)]';
if size(h0,1)==1
    Ng     =size(xg,1);
    yg     =y(round(Ny/2))*ones(Ng,1);
else
    xg     = [xg;xg];
    Ng     = size(xg,1);
    yg     =[ y(round(Ny/2))           *ones(Ng/2,1);...
             (y(round(Ny/2))-2*info.lc)*ones(Ng/2,1)];
end    
gauges = [xg,yg];
fid    = fopen(GAUGE,'w');
for jj = 1:size(gauges,1)
    fprintf(fid,'%i %i\n',gauges(jj,:));
end
info.gaugeFile = GAUGE;
info.Ng        = Ng;
info.gauges    = gauges;
fprintf(' info.gaugeFile = %s \n', info.gaugeFile)
fprintf(' info.Ng        = %u \n', info.Ng)
fprintf(' info.gauges    = %u %u \n', info.gauges(1,:))
fprintf('                   ...  \n', info.gauges(1,:))
%
save(info.fileName,'-struct','info')