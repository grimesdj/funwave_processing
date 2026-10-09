%% need variables:
% 0) requires the input bathymetry name as top-dir
jobName  = 'spreadRip';
runBATHYlist = {jobName};
% 1) need to know which cross-section to compare: {bathy, waves}
compare  = 'waves';% 'bathy' or 'waves'
varComp  = {'s'};
waveVars = {'h','t','s','d'};
%
% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
figDIR = [pwd,filesep,'../figures/'];
if ~exist(figDIR,'dir'), mkdir(pwd,'../figures'), end
%
% compile a list of all desired run_dirs
run_dirs = {};
for ii=1:length(runBATHYlist)
runBATHY = runBATHYlist{ii};
%
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
matDIR   = [runDIR,filesep,'mat_data'];
%
% the list of run directories are saved in:
tmp=load([matDIR,filesep,'runs_to_process.mat']);
run_dirs = cat(1,run_dirs,tmp.run_dirs);
end
clear tmp
%
%
% 
%% 2) now must deside to split based on bathy or waves:
grids = split(run_dirs,'_');
waves = grids(:,2);
grids = grids(:,1);
%
wave_str  = split(waves,waveVars);
if ~isempty(varComp)
    [~,idx]   = ismember(varComp,waveVars);
end
%
switch compare
  case 'bathy'
    run_list{1}  = run_dirs;
    %
    if ~isempty(varComp)
        tmp  = [grids(:,1), repmat({['_']},size(grids,1),1), ];
        tmp1 = [];
        for jj=1:length(varComp)
            tmp1 = cat(2,tmp1,repmat(varComp(jj),size(grids,1),1), wave_str(:,idx(jj)+1));
        end
        tmp = [tmp, tmp1];
    else
        tmp = [grids(:,1)];
    end
    run_names{1} = cellstr(cell2mat(tmp));
    run_label{1} = jobName;
    %
  case 'waves'
    %% 0) get a unique list of grids:
    ugrid = unique(grids);
    % 2) loop over ugrid
    for ii = 1:length(ugrid);
        % 3) find all runs with current ugrid
        idx1 = ismember(grids,ugrid(ii));
        run_list{ii} = run_dirs(idx1);
        %
        tmp = [];
        for jj=1:length(varComp)
            tmp = cat(2,tmp,repmat(varComp(jj),size(grids,1),1), wave_str(:,idx(jj)+1));
        end
        run_names{ii} = cellstr(cell2mat(tmp));
        run_label{ii} = ugrid{ii};
    end
end

%% 3) now loop over each group of runs... 
% plot:
% i)  Velocity/Vorticity: (Uex, Uex_avg, EKE, MKE) and (rms(VORT), rms(VORT_avg))
% ii) Momentum: (along/cross) shore at (Inner, Mid, Outer) surf
% iii) PGX vs ADX and RSX, and PGY vs ADY and RSX
% iv) diagnostic: PGY vs peak Uex
N    = length(run_list);
M    = length(run_list{1});
cm   = cmocean('thermal',M+1);
cm   = cm(1:M,:);
clrs = 1:M;
%
% figure (i) parameters:
% figure properties
xm = 2.5;
ym = 2.5;
pw = 4.5;
ph = 1.75;
ag = 0.25;
ppos1 = [xm ym           pw ph];
ppos2 = [xm ym+ph+ag     pw ph];
ppos3 = [xm ym+2*(ph+ag) pw ph];
cbpos0= [xm+pw+ag, ym, ag, ph];
ps    = [3*xm+3*ag+pw 2*ym+3*ag+3*ph];
%
fig0 = figure('units','centimeters');
fig0.Position(3:4)=ps;
fig0.PaperSize=ps;
fig0.PaperPosition=[0 0 ps];
% figure (ii) parameters:
xm = 2;
ym = 2;
pw = 3.5;
ph = 1.5;
ag = 0.2;
%
ppos31 = [xm           ym           pw ph];
ppos32 = [xm+pw+ag     ym           pw ph];
ppos33 = [xm+2*(pw+ag) ym           pw ph];
ppos21 = [xm           ym+ph+ag     pw ph];
ppos22 = [xm+pw+ag     ym+ph+ag     pw ph];
ppos23 = [xm+2*(pw+ag) ym+ph+ag     pw ph];
ppos11 = [xm           ym+2*(ph+ag) pw ph];
ppos12 = [xm+pw+ag     ym+2*(ph+ag) pw ph];
ppos13 = [xm+2*(pw+ag) ym+2*(ph+ag) pw ph];
cbpos1 = [xm+3*(pw+ag)+xm ym 2*ag 1.5*ph];
ps     = [2*xm+3*pw+8*ag 2*ym+3*ph+3*ag];
%
fig1 = figure('units','centimeters');
pos = get(fig1,'Position');
set(fig1,'Position',[pos(1:2) ps], 'Papersize',ps,'PaperPosition',[0 0 ps])
%
fig2 = figure('units','centimeters');
pos = get(fig2,'Position');
set(fig2,'Position',[pos(1:2) ps], 'Papersize',ps,'PaperPosition',[0 0 ps])
%
% figure (iii) parameters:
xm = 2;
ym = 2;
pw = 3;
ph = 3;
ag = 0.2;
%
f3ppos31 = [xm           ym           pw ph];
f3ppos32 = [xm+pw+ag     ym           pw ph];
f3ppos21 = [xm           ym+ph+ag     pw ph];
f3ppos22 = [xm+pw+ag     ym+ph+ag     pw ph];
f3ppos11 = [xm           ym+2*(ph+ag) pw ph];
f3ppos12 = [xm+pw+ag     ym+2*(ph+ag) pw ph];
cbpos3 = [xm+2*(pw+ag)+xm ym 2*ag 1.5*ph];
ps     = [2*xm+2*pw+6*ag 2*ym+3*ph+3*ag];
%
fig3 = figure('units','centimeters');
pos = get(fig3,'Position');
set(fig3,'Position',[pos(1:2) ps], 'Papersize',ps,'PaperPosition',[0 0 ps])
%
fig4 = figure('units','centimeters');
pos = get(fig4,'Position');
set(fig4,'Position',[pos(1:2) ps], 'Papersize',ps,'PaperPosition',[0 0 ps])
%
for ii=1:N
    run_dirs = run_list{ii};
    fprintf('\n plotting runs: \n')
    fprintf('\t\t %s \n', string(run_dirs))
    clf(fig0)
    figure(fig0)
    f0a1 = axes('units','centimeters','position',ppos1);
    f0a2 = axes('units','centimeters','position',ppos2);
    f0a3 = axes('units','centimeters','position',ppos3);    
    f0cb = axes('units','centimeters','position',cbpos0);
    clf(fig1)
    figure(fig1)
    for kk=1:3,
        for ll=1:3
            eval(['f1a',num2str(kk),num2str(ll),'=axes(''units'',''centimeters'',''position'',ppos',num2str(kk),num2str(ll),');'])
        end
    end
    f1cb = axes('units','centimeters','position',cbpos1);
    %
    clf(fig2)
    figure(fig2)
    for kk=1:3,
        for ll=1:3
            eval(['f2a',num2str(kk),num2str(ll),'=axes(''units'',''centimeters'',''position'',ppos',num2str(kk),num2str(ll),');'])
        end
    end
    f2cb = axes('units','centimeters','position',cbpos1);    
    %
    clf(fig3)
    figure(fig3)
    for kk=1:3,
        for ll=1:2
            eval(['f3a',num2str(kk),num2str(ll),'=axes(''units'',''centimeters'',''position'',f3ppos',num2str(kk),num2str(ll),');'])
        end
    end
    f3cb = axes('units','centimeters','position',cbpos3);    
    %
    clf(fig4)
    figure(fig4)
    for kk=1:3,
        for ll=1:2
            eval(['f4a',num2str(kk),num2str(ll),'=axes(''units'',''centimeters'',''position'',f3ppos',num2str(kk),num2str(ll),');'])
        end
    end
    f4cb = axes('units','centimeters','position',cbpos3);    
    %
    %
    for jj=1:M
        runID    = run_dirs{jj};
        infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
        info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
        %
        %
        fin = [info.rootMat,info.rootName,'dep.nc'];
        h = ncread(fin,'dep',[info.subDomain([1 3])] , [info.subDomain([2 4])]);
        %
        %% Velocity/Vorticity Stats:
        % estimate the energy/exchange for the mean fields:
        momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
        Umean   = ncread(momFile,'umean');
        Vmean   = ncread(momFile,'vmean');
        ETAmean = ncread(momFile,'etamean');
        Hmean   = max(h+ETAmean,0.1);
        %
        Tmean   = mean(Umean.*Hmean,[1 3],'omitnan');
        Us      = Tmean./mean(Hmean,[1 3],'omitnan');
        Umean   = mean(Umean-Us,3,'omitnan');
        Vmean   = mean(Vmean,3,'omitnan');
        ETAmean = mean(ETAmean,3,'omitnan');
        Hmean   = max(h+ETAmean,0.1);        
        %
        MKE     = sum(sqrt( Umean.^2 + Vmean.^2).*Hmean, 1,'omitnan')./sum(Hmean, 1,'omitnan');
        %
        tmp = Umean; tmp(Umean<0)=nan;
        tmp1 = Hmean; tmp1(Umean<0)=nan;
        Uex_avg = sum(tmp.*tmp1, 1, 'omitnan')./sum(tmp1, 1, 'omitnan');
        %
        [~,dUdy] = gradient(Umean./info.dy);
        [dVdx,~] = gradient(Vmean./info.dx);
        VORTavg  = dVdx-dUdy;
        %
        % read fields and calculate stats:
        x        = ncread(info.rotVelFile,'x');
        y        = ncread(info.rotVelFile,'y');
        t        = ncread(info.rotVelFile,'t');
        Urot     = ncread(info.rotVelFile,'Urot');
        Vrot     = ncread(info.rotVelFile,'Vrot');
        ETA      = ncread(info.rotVelFile,'eta');
        VORT     = ncread(info.rotVelFile,'VORT');
        disp('setting waterlevel to time-mean... bug in source code')
% $$$         H        = max(h+ETA,0.1);
        H = repmat(Hmean,1,1,size(Urot,3));
        %
        %
        % estimate energy/exchange statistics
        tmp      = sqrt( Urot.^2 + Vrot.^2 );
        EKE      = sum( tmp.*H , [1 3],'omitnan')./sum( H, [1 3],'omitnan');
        tmp      = Urot; tmp(Urot<0)=nan;
        tmp1     = H; tmp1(Urot<0)=nan;
        Uex      = sum( tmp.*tmp1, [1 3],'omitnan')./sum(tmp1, [1 3],'omitnan'); clear tmp tmp1
        %
        %
        figure(fig0), axes(f0a1), hold on,
        plot(x, rms(VORT,[1 3],'omitnan'),'-',x, rms(VORTavg,1,'omitnan'),':','color',cm(jj,:), 'linewidth',1)
        if jj==1,
            hl1 = legend({'$\mathrm{rms}(\bar{\omega})$','$\mathrm{rms}\langle\omega\rangle$'},'interpreter','latex','autoupdate','off','fontsize',8);
            hl1.AutoUpdate='off';
            hl1.ItemTokenSize=[10 10];
        end
        %
        axes(f0a2), hold on,
        plot(x, Uex,'-',x,Uex_avg,':','color',cm(jj,:), 'linewidth',1)
        if jj==1,
            hl2 = legend({'$U_\mathrm{ex}$','$\langle U\rangle_\mathrm{ex}$'},...
                         'interpreter','latex','autoupdate','off','fontsize',8);
            hl2.AutoUpdate='off';
            hl2.ItemTokenSize=[10 10];
        end
        %
        axes(f0a3), hold on,
        plot(x,EKE,'-',x,MKE,':','color',cm(jj,:), 'linewidth',1)
        if jj==1,
            hl3 = legend({'$U_\mathrm{eke}$','$U_\mathrm{mke}$'},...
                         'interpreter','latex','autoupdate','off','fontsize',8);
            hl3.AutoUpdate='off';
            hl3.ItemTokenSize=[10 10];
        end
        %
        %% Momentum Stats:
        momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
        % surfzone width and transect locations
% $$$         xSL = mean(info.x_shoreline);
% $$$         xBP = info.x_breakpoint;        
% $$$         iBP = find(info.x>=xBP,1,'first');
% $$$         Wsz = xBP-xSL;
% $$$         iINN= find(info.x>(xSL + Wsz/3),1,'first');
% $$$         iMID= find(info.x>(xSL + Wsz*2/3),1,'first');
% $$$         iOUT= find(info.x>(xSL + Wsz),1,'first');
        y0      = info.Ly/2;
        iOUT = find(info.x>=info.xc,1,'first');
        iINN = find(info.x>=50+(info.xc-50)/3,1,'first');
        iMID = find(info.x>=50+(info.xc-50)*2/3,1,'first');
        ylims   = [-500 500];
        yticks  = [-300 0 300];                
        %
        % 4.2.1) estimate time-averages of:
        %        cross-shore terms:
        %             advection, 
        tmp1 = ncread(momFile,'DxUUH');
        tmp2 = ncread(momFile,'DyUVH');
        ADX  = mean(tmp1,3,'omitnan')+mean(tmp2,3,'omitnan'); clear tmp1 tmp2
        %             pressure grad,
        PGX  = ncread(momFile,'PgrdX');
        PGX  = mean(PGX,3,'omitnan');
        %             radiation stress+BrkDissX,
        tmp1 = ncread(momFile,'DxSxx');        
        tmp2 = ncread(momFile,'DySxy');
        tmp3 = ncread(momFile,'BrkDissX');        
        RSX  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan') - mean(tmp3,3,'omitnan');
        clear tmp1 tmp2 tmp3
        %        along-shore terms:
        %             advection, 
        tmp1 = ncread(momFile,'DyVVH');
        tmp2 = ncread(momFile,'DxUVH');
        ADY  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan'); clear tmp1 tmp2
        %             pressure grad,
        PGY  = ncread(momFile,'PgrdY');
        PGY  = mean(PGY,3,'omitnan');
        %             radiation,
        tmp1 = ncread(momFile,'DySyy');
        tmp2 = ncread(momFile,'DxSxy');
        tmp3 = ncread(momFile,'BrkDissY');
        RSY  = mean(tmp1,3,'omitnan') + mean(tmp2,3,'omitnan') - mean(tmp3,3,'omitnan');
        clear tmp1 tmp2 tmp3
        %
        % mometnum terms need to be alongshore smoothed over ~25m (half-width of channel)
        Nflt = 25/info.dy;
        flt  = hamming(Nflt); flt = flt/sum(flt);
        ADX  = conv2(ADX,flt,'same');
        RSX  = conv2(RSX,flt,'same');
        PGX  = conv2(PGX,flt,'same');        
        %
        ADY  = conv2(ADY,flt,'same');
        RSY  = conv2(RSY,flt,'same');
        PGY  = conv2(PGY,flt,'same');        
        %
        %% Cross-shore Momentum terms:
        figure(fig1), 
        % plot the [inner, middle, outer]--Pressure Gradient
        axes(f1a11), hold on,
        plot(y-y0,PGX(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a12), hold on,
        plot(y-y0,PGX(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a13), hold on,
        plot(y-y0,PGX(:,iOUT),'-','color',cm(jj,:),'linewidth',1)        
        % plot the [inner, middle, outer]--Radiation Stress + Wave Breaking        
        axes(f1a21), hold on,
        plot(y-y0,RSX(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a22), hold on,
        plot(y-y0,RSX(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a23), hold on,
        plot(y-y0,RSX(:,iOUT),'-','color',cm(jj,:),'linewidth',1)
        % plot the [inner, middle, outer]--Advection
        axes(f1a31), hold on,
        plot(y-y0,ADX(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a32), hold on,
        plot(y-y0,ADX(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f1a33), hold on,
        plot(y-y0,ADX(:,iOUT),'-','color',cm(jj,:),'linewidth',1)
        %
        %
        %% Alongshore Momentum terms:
        figure(fig2), 
        % plot the [inner, middle, outer]--Pressure Gradient
        axes(f2a11), hold on,
        plot(y-y0,PGY(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a12), hold on,
        plot(y-y0,PGY(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a13), hold on,
        plot(y-y0,PGY(:,iOUT),'-','color',cm(jj,:),'linewidth',1)        
        % plot the [inner, middle, outer]--Radiation Stress + Wave Breaking        
        axes(f2a21), hold on,
        plot(y-y0,RSY(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a22), hold on,
        plot(y-y0,RSY(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a23), hold on,
        plot(y-y0,RSY(:,iOUT),'-','color',cm(jj,:),'linewidth',1)
        % plot the [inner, middle, outer]--Advection
        axes(f2a31), hold on,
        plot(y-y0,ADY(:,iINN),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a32), hold on,
        plot(y-y0,ADY(:,iMID),'-','color',cm(jj,:),'linewidth',1)
        axes(f2a33), hold on,
        plot(y-y0,ADY(:,iOUT),'-','color',cm(jj,:),'linewidth',1)
        %
        %% Cross-shore Momentum: PGX vs RSX & ADX
        figure(fig3), 
        % plot the [inner]--
        axes(f3a31), hold on,
        plot(RSX(:,iINN), PGX(:,iINN),'.','color',cm(jj,:),'markersize',1)
        axes(f3a32), hold on,
        plot(ADX(:,iINN),PGX(:,iINN),'.','color',cm(jj,:),'markersize',1)
        % plot the [middle]--
        axes(f3a21), hold on,
        plot(RSX(:,iMID), PGX(:,iMID),'.','color',cm(jj,:),'markersize',1)
        axes(f3a22), hold on,
        plot(ADX(:,iMID), PGX(:,iMID),'.','color',cm(jj,:),'markersize',1)
        % plot the [outer]--
        axes(f3a11), hold on,
        plot(RSX(:,iINN), PGX(:,iMID),'.','color',cm(jj,:),'markersize',1)
        axes(f3a12), hold on,
        plot(ADX(:,iOUT), PGX(:,iOUT),'.','color',cm(jj,:),'markersize',1)
        %
        %
        %% Alongshore Momentum terms:
        figure(fig4),
        % plot the [inner]--
        axes(f4a31), hold on,
        plot(RSY(:,iINN), PGY(:,iINN),'.','color',cm(jj,:),'markersize',1)
        axes(f4a32), hold on,
        plot(ADY(:,iINN),PGY(:,iINN),'.','color',cm(jj,:),'markersize',1)
        % plot the [middle]--
        axes(f4a21), hold on,
        plot(RSY(:,iMID), PGY(:,iMID),'.','color',cm(jj,:),'markersize',1)
        axes(f4a22), hold on,
        plot(ADY(:,iMID), PGY(:,iMID),'.','color',cm(jj,:),'markersize',1)
        % plot the [outer]--
        axes(f4a11), hold on,
        plot(RSY(:,iINN), PGY(:,iMID),'.','color',cm(jj,:),'markersize',1)
        axes(f4a12), hold on,
        plot(ADY(:,iOUT), PGY(:,iOUT),'.','color',cm(jj,:),'markersize',1)
        %
        %
    end
    %% Velocity/Vorticity Statistics
    figure(fig0)
    axes(f0cb)
    imagesc(0,clrs,reshape(cm,M,1,3))
    set(f0cb,'ydir','normal','yaxislocation','right','xaxislocation','top','fontsize',6,'ytick',1:M,'yticklabel',run_names{ii},'ticklabelinterpreter','latex','xtick',[])
    xlabel(f0a1,'$x$ [m]','interpreter','latex')
    ylabel(f0a1,'(s$^{-1}$)','interpreter','latex')
    set(f0a1,'tickdir','out','ticklabelinterpreter','latex')
    ylabel(f0a2,'(m/s)','interpreter','latex')
    set(f0a2,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    ylabel(f0a3,'(m/s)','interpreter','latex')
    set(f0a3,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    if ~isempty(varComp)
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_param_',char(varComp),'_velocity_and_vorticity_stats.pdf'];
    else
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_velocity_and_vorticity_stats.pdf'];        
    end
    exportgraphics(fig0,figname)
    %
    %% Cross-shore Momentum
    figure(fig1)
    axes(f1cb)
    imagesc(0,clrs,reshape(cm,M,1,3))
    set(f1cb,'ydir','normal','yaxislocation','right','xaxislocation','top','fontsize',6,'ytick',1:M,'yticklabel',run_names{ii},'ticklabelinterpreter','latex','xtick',[])
    xlabel(f1a32,'$y$ [m]','interpreter','latex')
    ylabel(f1a21,'[m/s]$^2 \times 10^{-2}$','interpreter','latex')
    title(f1a11,'Inner Surfzone','interpreter','latex','fontsize',8)
    title(f1a12,'Middle Surfzone','interpreter','latex','fontsize',8)
    title(f1a13,'Outer Surfzone','interpreter','latex','fontsize',8)
    ylabel(f1a13,'Pressure','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f1a23,'Wave','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f1a33,'Advection','rotation',0,'horizontalalignment','left','fontsize',8)
    set(f1a31,'tickdir','out','ticklabelinterpreter','latex')
    set([f1a21 f1a11],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    set([f1a32 f1a33],'tickdir','out','ticklabelinterpreter','latex','yticklabel',[])
    set([f1a12 f1a13 f1a22 f1a23],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
    set([f1a13 f1a23 f1a33],'yaxislocation','right')
    for kk=1:3
        for ll=1:3
            eval(['f1a',num2str(kk),num2str(ll),'.YAxis.Exponent=-2;'])
            eval(['f1a',num2str(kk),num2str(ll),'.YAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f1a',num2str(kk),num2str(ll),'.YAxis.TickValues=[-1 0 1]*1e-2;'])                                    
            eval(['f1a',num2str(kk),num2str(ll),'.XAxis.Limits=[',num2str(ylims),'];'])
            eval(['f1a',num2str(kk),num2str(ll),'.XAxis.TickValues=[',num2str(yticks),'];'])            
            eval(['f1a',num2str(kk),num2str(ll),'.FontSize=8;'])
            eval(['f1a',num2str(kk),num2str(ll),'.Box=''on'';'])
            eval(['grid(f1a',num2str(kk),num2str(ll),',''on'');'])                                                                                    
        end
    end
    if ~isempty(varComp)
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_param_',char(varComp),'_dominant_x_momentum.pdf'];
    else
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_dominant_x_momentum.pdf'];        
    end
    exportgraphics(fig1,figname)
    %
    %% Alongshore Momentum
    figure(fig2)
    axes(f2cb)
    imagesc(0,clrs,reshape(cm,M,1,3))
    set(f2cb,'ydir','normal','yaxislocation','right','xaxislocation','top','fontsize',6,'ytick',1:M,'yticklabel',run_names{ii},'ticklabelinterpreter','latex','xtick',[])
    xlabel(f2a32,'$y$ [m]','interpreter','latex')
    ylabel(f2a21,'[m/s]$^2 \times 10^{-2}$','interpreter','latex')
    title(f2a11,'Inner Surfzone','interpreter','latex','fontsize',8)
    title(f2a12,'Middle Surfzone','interpreter','latex','fontsize',8)
    title(f2a13,'Outer Surfzone','interpreter','latex','fontsize',8)
    ylabel(f2a13,'Pressure','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f2a23,'Wave','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f2a33,'Advection','rotation',0,'horizontalalignment','left','fontsize',8)
    set(f2a31,'tickdir','out','ticklabelinterpreter','latex')
    set([f2a21 f2a11],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    set([f2a32 f2a33],'tickdir','out','ticklabelinterpreter','latex','yticklabel',[])
    set([f2a12 f2a13 f2a22 f2a23],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
    set([f2a13 f2a23 f2a33],'yaxislocation','right')
    for kk=1:3
        for ll=1:3
            eval(['f2a',num2str(kk),num2str(ll),'.YAxis.Exponent=-2;'])
            eval(['f2a',num2str(kk),num2str(ll),'.YAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f2a',num2str(kk),num2str(ll),'.YAxis.TickValues=[-1 0 1]*1e-2;'])                                    
            eval(['f2a',num2str(kk),num2str(ll),'.XAxis.Limits=[',num2str(ylims),'];'])
            eval(['f2a',num2str(kk),num2str(ll),'.XAxis.TickValues=[',num2str(yticks),'];'])                        
            eval(['f2a',num2str(kk),num2str(ll),'.FontSize=8;'])
            eval(['f2a',num2str(kk),num2str(ll),'.Box=''on'';'])
            eval(['grid(f2a',num2str(kk),num2str(ll),',''on'');'])                                                                        
        end
    end
    if ~isempty(varComp)
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_param_',char(varComp),'_dominant_y_momentum.pdf'];
    else
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_dominant_y_momentum.pdf'];        
    end
    exportgraphics(fig2,figname)
    %
    %% Cross-shore Momentum: PGX vs ...
    figure(fig3)
    axes(f3cb)
    imagesc(0,clrs,reshape(cm,M,1,3))
    set(f3cb,'ydir','normal','yaxislocation','right','xaxislocation','top','fontsize',6,'ytick',1:M,'yticklabel',run_names{ii},'ticklabelinterpreter','latex','xtick',[])
    xlabel(f3a31,'Wave [m/s]$^2 \times 10^{-2}$','interpreter','latex')
    xlabel(f3a32,'Advection','interpreter','latex')    
    ylabel(f3a21,'Pressure  [m/s]$^2 \times 10^{-2}$','interpreter','latex')
% $$$     title(f1a11,'Inner Surfzone','interpreter','latex','fontsize',8)
% $$$     title(f1a12,'Middle Surfzone','interpreter','latex','fontsize',8)
% $$$     title(f1a13,'Outer Surfzone','interpreter','latex','fontsize',8)
    ylabel(f3a12,'Outer-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f3a22,'Mid-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f3a32,'Inner-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    set(f3a31,'tickdir','out','ticklabelinterpreter','latex')
    set([f3a21 f3a11],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    set( f3a32,'tickdir','out','ticklabelinterpreter','latex','yticklabel',[])
    set([f3a12 f3a22],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
    set([f3a12 f3a22 f3a32],'yaxislocation','right')
    for kk=1:3
        for ll=1:2
            eval(['f3a',num2str(kk),num2str(ll),'.YAxis.Exponent=-2;'])
            eval(['f3a',num2str(kk),num2str(ll),'.YAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f3a',num2str(kk),num2str(ll),'.YAxis.TickValues=[-1 0 1]*1e-2;'])
            eval(['f3a',num2str(kk),num2str(ll),'.XAxis.Exponent=-2;'])                        
            eval(['f3a',num2str(kk),num2str(ll),'.XAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f3a',num2str(kk),num2str(ll),'.XAxis.TickValues=[-1 1]*1e-2;'])            
            eval(['f3a',num2str(kk),num2str(ll),'.FontSize=8;'])
            eval(['f3a',num2str(kk),num2str(ll),'.Box=''on'';'])
            eval(['grid(f3a',num2str(kk),num2str(ll),',''on'');'])                                                                                    
        end
    end
    if ~isempty(varComp)
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_param_',char(varComp),'_PGX_vs_RSX_ADX.pdf'];
    else
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_PGX_vs_RSX_ADX.pdf'];        
    end
    exportgraphics(fig3,figname)
    %
    %% Alongshore Momentum: PGY vs ...
    figure(fig4)
    axes(f4cb)
    imagesc(0,clrs,reshape(cm,M,1,3))
    set(f4cb,'ydir','normal','yaxislocation','right','xaxislocation','top','fontsize',6,'ytick',1:M,'yticklabel',run_names{ii},'ticklabelinterpreter','latex','xtick',[])
    xlabel(f4a31,'Wave [m/s]$^2 \times 10^{-2}$','interpreter','latex')
    xlabel(f4a32,'Advection','interpreter','latex')    
    ylabel(f4a21,'Pressure  [m/s]$^2 \times 10^{-2}$','interpreter','latex')
% $$$     title(f1a11,'Inner Surfzone','interpreter','latex','fontsize',8)
% $$$     title(f1a12,'Middle Surfzone','interpreter','latex','fontsize',8)
% $$$     title(f1a13,'Outer Surfzone','interpreter','latex','fontsize',8)
    ylabel(f4a12,'Outer-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f4a22,'Mid-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    ylabel(f4a32,'Inner-Surfzone','rotation',0,'horizontalalignment','left','fontsize',8)
    set(f4a31,'tickdir','out','ticklabelinterpreter','latex')
    set([f4a21 f4a11],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[])
    set( f4a32,'tickdir','out','ticklabelinterpreter','latex','yticklabel',[])
    set([f4a12 f4a22],'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],'yticklabel',[])
    set([f4a12 f4a22 f4a32],'yaxislocation','right')
    for kk=1:3
        for ll=1:2
            eval(['f4a',num2str(kk),num2str(ll),'.YAxis.Exponent=-2;'])
            eval(['f4a',num2str(kk),num2str(ll),'.YAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f4a',num2str(kk),num2str(ll),'.YAxis.TickValues=[-1 0 1]*1e-2;'])
            eval(['f4a',num2str(kk),num2str(ll),'.XAxis.Exponent=-2;'])            
            eval(['f4a',num2str(kk),num2str(ll),'.XAxis.Limits=[-1.5 1.5]*1e-2;'])
            eval(['f4a',num2str(kk),num2str(ll),'.XAxis.TickValues=[-1 0 1]*1e-2;'])            
            eval(['f4a',num2str(kk),num2str(ll),'.FontSize=8;'])
            eval(['f4a',num2str(kk),num2str(ll),'.Box=''on'';'])
            eval(['grid(f4a',num2str(kk),num2str(ll),',''on'');'])                                                                                    
        end
    end
    if ~isempty(varComp)
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_param_',char(varComp),'_PGY_vs_RSY_ADY.pdf'];
    else
        figname = [figDIR,'compare_',run_label{ii},'_',compare,'_PGY_vs_RSY_ADY.pdf'];        
    end
    exportgraphics(fig4,figname)
    %
end
