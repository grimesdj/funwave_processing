function make_1panel_video(vidName,x,y,t,data,alims,clims,clr_map,label1,label2)

%%    ii) make a movie
% define limits 
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean(clr_map);
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 3;
open(vid)
%
Nt = length(t);
for jj = 1:Nt;
    % plot avg
    imagesc(ax0,y,x,squeeze(data(:,:,jj)')), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    xlabel(ax0,'$y$ [m]','interpreter','latex')
    ylabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf(['$t$ = %1.1f min, ',label1],(mean(t(jj))-t(1))/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(1:2),'xlim',alims(3:4)+y(1))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,label2,'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    if vid.FrameCount==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close(vid)
close(fig)
