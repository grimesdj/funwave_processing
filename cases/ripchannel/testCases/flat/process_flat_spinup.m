clear all
close all

addpath(genpath('~/git/funwave/code/'))
%% compile array of eta during spin-up for flat bottom run. Trying to understand interaction with the sponges.

%% First create the "info" structure that has directory and file information for code:
info = struct([]);
info(1).runName   = 'spinup';
info.rootName  = info.runName;
info.gaugeFile = '/scratch/grimesdj/ripchannel/flat/gauge.txt';
info.bathyFile = '/scratch/grimesdj/ripchannel/flat/flatBathy.mat';
info.rootSim   = '/scratch/grimesdj/ripchannel/flat/spinup/';
info.rootMat   = '/scratch/grimesdj/ripchannel/flat/mat_data_spinup/';
info.fileName  = [info.rootMat,filesep,info.rootName,'_info.mat'];

if ~exist(info.rootMat,'dir')
    eval(['!mkdir -p ',info.rootMat])
end

%% Read input spectrum from "LOG.txt"
fin = [info.rootSim,filesep,'LOG.txt'];
[freq, dire, amp] = parse_wavemaker_info_from_LOG(fin,0);
if isempty(freq)
    [freq, dire, amp] = parse_wavemaker_info_LOG_aggregated(fin,0);
    if isempty(freq)
        disp('no wavemaker informamtion read from LOG.txt')
    end
end
                 
info.freq=freq;
info.dire=dire;
info.amp =amp;

%% Read final Hs field (Mglob and Nglob from input file)
Mglob = 1364;
Nglob = 2999;
HsFiles = dir([info.rootSim,filesep,'output',filesep,'Hsig_*']);
info.Hs = 0;
Nf      = length(HsFiles);
for jj=1:Nf
    fin = [HsFiles(jj).folder,filesep,HsFiles(jj).name];
    fid = fopen(fin);
    dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
    info.Hs = info.Hs+dum';% so I'm transposing so rows are y-coord and columns are x-coord
    fclose(fid);
end
info.Hs = info.Hs/Nf;

etaFiles = dir([info.rootSim,filesep,'output',filesep,'eta_*']);
Nf       = length(etaFiles);
eta = [];
%
%%    ii) make a movie
x = [0:(Mglob-1)]*0.5;
y = [0:(Nglob-1)]*1.0;
% define limits 
alims = [x(1) x(end) y(1) y(end)];
clims = [-1 1]*1e-1;
clrs  = clims(1):diff(clims)/255:clims(2);
cm    = cmocean('balance');
[fig,ax0,ax00,cx01,ps,ppos,pos] = get_1panel_video_figure_info(alims);
delete(ax00)
%
%
vidName= [info.rootMat,info.runName,'_eta'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
vid.Quality  = 100;
vid.FrameRate= 5;
open(vid)
%
clims = [-0.5 0.5]
clrs  = clims(1):diff(clims)/255:clims(2);
cm = cmocean('balance');
%
for jj=1:Nf
    fin = [etaFiles(jj).folder,filesep,etaFiles(jj).name];
    fid = fopen(fin);
    dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
    eta = cat(3,eta,dum');
    fclose(fid);
    %
    imagesc(ax0,x,y,dum'), 
    caxis(ax0,clims)
    colormap(ax0,cm)
    ylabel(ax0,'$y$ [m]','interpreter','latex')
    xlabel(ax0,'$x$ [m]','interpreter','latex')
    title_str = sprintf('$t$ = %1.1f min, $\\langle \\omega \\rangle$ ',jj/60);
    set(ax0,'tickdir','out','ticklabelinterpreter','latex','fontsize',25,'ydir','normal','color',0.8*[1 1 1],'xdir','reverse','ylim',alims(3:4),'xlim',alims(1:2))
    title(ax0,title_str,'interpreter','latex','fontsize',15,'horizontalalignment','left','units','normalized','position',[0.01 1.1 0]) 
    %
    % make colorbar
    imagesc(cx01,clrs,0,reshape(cm,1,256,3))
    xlabel(cx01,{'(s$^{-1}$)'},'interpreter','latex','rotation',0)%,'horizontalalignment','right')
    set(cx01,'ytick',[],'xaxislocation','top','tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    %
    % 
    % get+write frame
    set(fig,'position',pos);
    frame = getframe(fig);
    try FC = vid.FrameCount;
    catch
        FC=0;
    end
    if FC==0
        fsize = size(frame.cdata,1,2);
    elseif any(size(frame.cdata,1,2)~=fsize)
        frame.cdata = imresize(frame.cdata,fsize);
    end
    writeVideo(vid,frame)
    %
end
close all
info = process_FUNWAVE_virtual_moorings(info);

save('/scratch/grimesdj/ripchannel/flat/mat_data_spinup/flat_spinup_run.mat','-v7.3')
