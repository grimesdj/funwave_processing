%% create compendium plot of exchange and enstrophy for a set of runs:
clear all
close all
% 
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
%
%% 0) where are we looking? archiving?
rootDIR  = '/scratch/grimesdj/ripchannel/';
figDIR   = '/scratch/grimesdj/ripchannel/figures/';
outDIR   = '/scratch/grimesdj/ripchannel/mat_data/';
%
if ~exist(figDIR,'dir'), eval(['!mkdir ', figDIR]), end
if ~exist(outDIR,'dir'), eval(['!mkdir ', outDIR]), end
%
%% 1) what cases to process?
CASES = {'spreadRip','gapRip','gapRip1','highRip','highRip1'};
%
%% 1.01) plot as we go?
plotter = 1;% 1=make plots, 0=don't;
% initialize w/ zero 
plotFlag= 0;
axisNum = 0;
% for velocity plots (maximum of 5-panels)
xm = 2;
ym = 2;
pw = 2.5;% 400 in x
ph = 6;% +/- 500 in y
ag = 0.5;
%% 1.1) Preallocate archive variables:
% INFO = struct([]);
Uex = []; Uex_mean = []; Uex_eddy = [];
Uex_vs_x = []; Uex_mean_vs_x = []; Uex_eddy_vs_x = [];
Urms= []; Urms_mean= []; Urms_eddy= [];
ENS = []; ENS_mean = []; ENS_eddy = [];
dETA= []; Uscale   = [];
Umax= []; Umax_mean= []; Umax_eddy= [];% all at bar crest
Umax_BP= []; Umax_mean_BP= []; Umax_eddy_BP= [];% at avg break-point
% time-series of "maximum" speed at crest and break-point
Umax_vs_t      = []; Umax_vs_t_BP      = [];
Umax_eddy_vs_t = []; Umax_eddy_vs_t_BP = [];                
% mrkr= {}; sz = {}; cblbl = {}; cbttl = [];
Uex_channel_vs_x = []; Uex_ambient_vs_x = [];
Uex_channel = [];      Uex_ambient = [];
Uex_mean_channel_vs_x = []; Uex_mean_ambient_vs_x = [];
Uex_mean_channel = []; Uex_mean_ambient = [];
Uex_eddy_channel_vs_x = []; Uex_eddy_ambient_vs_x = [];
Uex_eddy_channel = []; Uex_eddy_ambient = [];
Urms_channel     = []; Urms_ambient = [];
Urms_mean_channel= []; Urms_mean_ambient = [];
Urms_eddy_channel= []; Urms_eddy_ambient = [];
ENS_channel      = []; ENS_ambient = [];
ENS_mean_channel = []; ENS_mean_ambient = [];
ENS_eddy_channel = []; ENS_eddy_ambient = [];
Vavg_trough      = []; Vrms_trough      = [];
Uavg_crest       = []; Urms_crest       = [];
Uavg_BP          = []; Urms_BP          = [];
height = []; period = []; spread = []; direction = [];
bar_location = []; bar_width = []; bar_amplitude = [];
channel_length = []; channel_amplitude_ratio = [];
%
%% 2) loop over the simulations,
iter = 0;
for cc = 1:length(CASES)
    CASE = CASES{cc};
    runDIR   = [rootDIR,CASE];
    matDIR   = [runDIR,filesep,'mat_data'];
    %
    % the list of run directories are saved in:
    load([matDIR,filesep,'runs_to_process.mat'])
    % brings in cell array: run_dirs
    % for example,
    % run_dirs =
    %   6x1 cell array
    %    {'planar1D_h05t08s00d00'}
    %    {'planar1D_h05t10s00d00'}
    %    {'planar1D_h10t08s00d00'}
    %    {'planar1D_h10t10s00d00'}
    %    {'planar1D_h15t08s00d00'}
    %    {'planar1D_h15t10s00d00'}
    Ndirs = length(run_dirs);
    %% 3) loop over each run in simulation,
    for rr = 1:length(run_dirs);
        iter = iter+1;
        RUN  = run_dirs{rr};
        %% 4) load info file:
        fprintf('\n processing: %s %s \n', CASE,RUN)    
        infoFile = dir([matDIR,filesep,'*','info','*',RUN,'.mat']);
        info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
        INFO(iter)=info;
        
        %% 4.1) get grid & wave parameters for run:
        grid_waves = split(info.runName,'_');
        grid  = grid_waves{1};
        waves = split(grid_waves{2},{'h','t','s','d'});
        %
        height(iter)    = str2num(waves{2})/10;
        period(iter)    = str2num(waves{3});
        spread(iter)    = str2num(waves{4});
        direction(iter) = str2num(waves{5});
        %
        bar_location(iter)            = info.xc;
        bar_width(iter)               = info.wc;
        bar_amplitude(iter)           = info.ac;
        channel_length(iter)          = info.lc;
        channel_amplitude_ratio(iter) = info.rc;

        %% create plot title string:
        title_string{1} = sprintf(' %s, $x_c=$%3.0f m, $l_c=$%3.0f m', grid(1:3), info.xc,info.lc);
        title_string{2,1} = sprintf(' $H=$%1.1f m, $T=$%2.0f s, $\\sigma_\\theta=$%d, $D=$%d ', height(iter),period(iter),spread(iter),direction(iter));

        %% 4.1.1) define label sizes/shapes/colors?
       
        %% 4.2) load post-processed variables:
        rotFile = [info.rootMat,'funwave_',RUN,'_velocity_decomposition.nc'];
        % get (x,y,dep)
        x   = ncread(rotFile,'x');
        y   = ncread(rotFile,'y');
        dep = ncread(rotFile,'dep');
        % get u=<u>+\bar{u};
        Umean = ncread(rotFile,'Urot_mean');
        Vmean = ncread(rotFile,'Vrot_mean');
        U     = ncread(rotFile,'Urot');
        V     = ncread(rotFile,'Vrot');
        % get <\eta>
        ETA   = ncread(rotFile,'eta');
        % mask WL below z=zbottom;
        ETA(ETA>dep) = 0;
        %
        % get vorticity
        VORT      = ncread(rotFile,'VORT');
        VORT_mean = curl(x,y,Umean,Vmean);
        
        %% 4.3) compute/compile statistics for run:
        % create depth mask (min-depth-resolved=0.01m, min-depth-normalize=0.1m)
        H         = dep+ETA;
        mask      = H>0.01;
        H(~mask)  = 0;
        Hmean     = max( mean(H,3), 0.1);
        %
        % mask of "offshore" directed flow for exchange calculations
        iOFF      = (Umean+U)>0;
        iOFF_mean = (Umean)>0;
        iOFF_eddy = (U)>0;
        %
        % !If domains are different, need to adjust array sizes after iter==1!
        Uex_vs_x(:,iter)       = mean( (Umean+U).*iOFF     .*H, [1 3]) ./ mean(Hmean, [1 3]);
        Uex_mean_vs_x(:,iter)  = mean( (Umean  ).*iOFF_mean.*H, [1 3]) ./ mean(Hmean, [1 3]);
        Uex_eddy_vs_x(:,iter)  = mean( (U      ).*iOFF_eddy.*H, [1 3]) ./ mean(Hmean, [1 3]);
        %
        Urms     (:,iter)      = sqrt(mean(sum( (Umean+U).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
        Urms_mean(:,iter)      = sqrt(mean(sum( (Umean  ).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
        Urms_eddy(:,iter)      = sqrt(mean(sum( (U      ).^2.*H, 1,'omitnan'), 3,'omitnan') ./ sum(Hmean, 1,'omitnan'));
        %
        ENS(:,iter)       = mean( (VORT_mean+VORT).^2.*H, [1 3])./mean( Hmean, 1);
        ENS_mean(:,iter)  = mean( (VORT_mean     ).^2.*H, [1 3])./mean( Hmean, 1);
        ENS_eddy(:,iter)  = mean( (VORT          ).^2.*H, [1 3])./mean( Hmean, 1);
        %
        %% 5) get locations for transects/bulk stats for each:
        iX0 = (x>0.5*(info.xc + mean(info.x_shoreline)) & x<info.xc)';
        iX1 = find(x>=info.xc,1,'first');
        iX2 = find(x>=info.x_breakpoint,1,'first');
        %
        iY0 =  y<info.Ly/2 - 5*info.lc | y>info.Ly/2 + 5*info.lc;
        iY1 = (y>info.Ly/2 -   info.lc & y<info.Ly/2 +   info.lc);

        %% 5.1) extract/archive statistics:
        Uex_channel_vs_x     (:,iter)  = mean( (Umean(iY1,:,:)+U(iY1,:,:)).*iOFF(iY1,:,:)     .*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
        Uex_mean_channel_vs_x(:,iter)  = mean( (Umean(iY1,:,:)           ).*iOFF_mean(iY1,:,:).*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
        Uex_eddy_channel_vs_x(:,iter)  = mean( (U(iY1,:,:)               ).*iOFF_eddy(iY1,:,:).*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]);
        %
        Uex_ambient_vs_x     (:,iter)  = mean( (Umean(iY0,:,:)+U(iY0,:,:)).*iOFF(iY0,:,:)     .*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
        Uex_mean_ambient_vs_x(:,iter)  = mean( (Umean(iY0,:,:)           ).*iOFF_mean(iY0,:,:).*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
        Uex_eddy_ambient_vs_x(:,iter)  = mean( (U(iY0,:,:)               ).*iOFF_eddy(iY0,:,:).*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]);
        %
        Uex(iter)       = Uex_vs_x(iX1,iter);
        Uex_mean(iter)  = Uex_mean_vs_x(iX1,iter);
        Uex_eddy(iter)  = Uex_eddy_vs_x(iX1,iter);
        %
        Uex_channel     (iter)  = Uex_channel_vs_x(iX1,iter);
        Uex_mean_channel(iter)  = Uex_mean_channel_vs_x(iX1,iter);
        Uex_eddy_channel(iter)  = Uex_eddy_channel_vs_x(iX1,iter);
        %
        Uex_ambient     (iter)  = Uex_ambient_vs_x(iX1,iter);
        Uex_mean_ambient(iter)  = Uex_mean_ambient_vs_x(iX1,iter);
        Uex_eddy_ambient(iter)  = Uex_eddy_ambient_vs_x(iX1,iter);
        %
        %
        Urms_channel     (:,iter)  = sqrt(mean( (Umean(iY1,:,:)+U(iY1,:,:)).^2.*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]));
        Urms_mean_channel(:,iter)  = sqrt(mean( (Umean(iY1,:,:)           ).^2.*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]));
        Urms_eddy_channel(:,iter)  = sqrt(mean( (U(iY1,:,:)               ).^2.*H(iY1,:,:), [1 3]) ./ mean(Hmean(iY1,:,:), [1 3]));
        %
        Urms_ambient     (:,iter)  = sqrt(mean( (Umean(iY0,:,:)+U(iY0,:,:)).^2.*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]));
        Urms_mean_ambient(:,iter)  = sqrt(mean( (Umean(iY0,:,:)           ).^2.*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]));
        Urms_eddy_ambient(:,iter)  = sqrt(mean( (U(iY0,:,:)               ).^2.*H(iY0,:,:), [1 3]) ./ mean(Hmean(iY0,:,:), [1 3]));
        %
        ENS_channel     (:,iter)  = mean( (VORT_mean(iY1,:,:)+VORT(iY1,:,:)).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
        ENS_mean_channel(:,iter)  = mean( (VORT_mean(iY1,:,:)              ).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
        ENS_eddy_channel(:,iter)  = mean( (VORT(iY1,:,:)                   ).^2.*H(iY1,:,:), [1 3])./mean( Hmean(iY1,:,:), 1);
        %
        ENS_ambient     (:,iter)  = mean( (VORT_mean(iY0,:,:)+VORT(iY0,:,:)).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
        ENS_mean_ambient(:,iter)  = mean( (VORT_mean(iY0,:,:)              ).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
        ENS_eddy_ambient(:,iter)  = mean( (VORT(iY0,:,:)                   ).^2.*H(iY0,:,:), [1 3])./mean( Hmean(iY0,:,:), 1);
        %
        deta= mean(ETA(iY1,:),1)-mean(ETA(iY0,:),1);
        dETA     (:,iter)= deta;
        Uscale   (iter)  = real( sqrt(-2*9.8*mean(deta(iX0))));
        Umax     (iter)  = max(U(iY1,iX1,:)+Umean(iY1,iX1,:), [], [1 3]);
        Umax_mean(iter)  = max(             Umean(iY1,iX1,:), [], [1 3]);
        Umax_eddy(iter)  = max(U(iY1,iX1,:)                 , [], [1 3]);        
        Umax_vs_t     (:,iter)  = max(U(iY1,iX1,:)+Umean(iY1,iX1,:), [], 1);
        Umax_eddy_vs_t(:,iter)  = max(U(iY1,iX1,:)                 , [], 1);
        %
        Umax_BP     (iter)  = max(U(iY1,iX2,:)+Umean(iY1,iX2,:), [], [1 3]);
        Umax_mean_BP(iter)  = max(             Umean(iY1,iX2,:), [], [1 3]);
        Umax_eddy_BP(iter)  = max(U(iY1,iX2,:)                 , [], [1 3]);
        Umax_vs_t_BP     (:,iter)  = max(U(iY1,iX2,:)+Umean(iY1,iX2,:), [], 1);
        Umax_eddy_vs_t_BP(:,iter)  = max(U(iY1,iX2,:)                 , [], 1);
        %
        Vmean_trough(:,iter) = mean( Vmean(:,iX0).*Hmean(:,iX0), 2 ) ./ mean(Hmean(:,iX0), 2);
        Veddy_trough(:,iter) = sqrt(mean( (V(:,iX0,:).^2).*H(:,iX0,:), [2 3]) ./ mean(Hmean(:,iX0), 2));
        %
        Umean_crest(:,iter) =  Umean(:,iX1);
        Ueddy_crest(:,iter) = sqrt(mean( (U(:,iX1,:).^2).*H(:,iX1,:), [2 3]) ./ Hmean(:,iX1));
        %
        Umean_BP(:,iter) =  Umean(:,iX2);
        Ueddy_BP(:,iter) = sqrt(mean( (U(:,iX2,:).^2).*H(:,iX2,:), [2 3]) ./ Hmean(:,iX2));

        
        %% 5.2) estimate Moulton scaling,
        %% Kludge to find bar-crest:
        icrest = find(x>=bar_location(iter),1,'first');
        hcrest = 0.03*(x(icrest)-50)-bar_amplitude(iter).*exp(-0.5.*((x(icrest)-bar_location(iter))./bar_width(iter)).^2);
        %-------------------------------------
        %% estimate H_breaking:
        %       Hbr^2 cg_br          = H0^2 cg0;
        %  Shallow water at breaking:
        %       Hbr^2 \sqrt(g h_br)  = H0^2 cg0;
        %  Assume: hbr-->Hbr/\gamma_br:
        %       Hbr^5 (g/\gamma_br)  = H0^4 cg0^2; 
        %       Hbr                  =(H0^4 cg0^2 \gamma_br/g)^(1/5)
        g = 9.8;
        gam = 0.35;
        dwdk0 = groupSpeed_FunwaveTVD(2*pi/period(iter),info.hWM);
        Hbr0 = nthroot( height(iter).^4.*dwdk0.^2.*(gam/g), 5);
        dETA_moulton0 = -gam.^2/16.*(cosd(direction(iter)).^2 + 0.5).*max(Hbr0/gam-hcrest,0);
        Uscale0(iter) = sqrt(-2*g*dETA_moulton0);
        
        %% 5.3) log anything else?

        %% 6) make diagnostic plots?
        if plotter
            if axisNum==0
                %% check how many runs have current grid:
                grids = split(run_dirs,'_');
                grids = char(grids(:,1));
                idx   = ismember(grids(:,1:end-1),grid(1:end-1),'rows');
                Nax   = sum(idx);
                %
                fprintf('\t preparing %d-axes for grids: %sX \n', Nax, grid(1:end-1))    
                for pp = 1:Nax
                    eval(sprintf('ppos%d = [xm+%d*(pw+ag) ym pw ph];',pp,pp-1));
                end
                ps    = [2*xm+Nax*(pw+ag)+6*ag 2*ym+ph];
                cbpos = [xm+Nax*(pw+ag)+ag, ym, ag, ph/2];
                %
                fig0  = figure('units','centimeters');
                fig00 = figure('units','centimeters');
                fig000 = figure('units','centimeters');
                fig0000 = figure('units','centimeters');
                fig00000 = figure('units','centimeters');
                fig0.Position(3:4)  = ps;
                fig00.Position(3:4) = ps;
                fig000.Position(3:4) = ps;
                fig0000.Position(3:4) = ps;
                fig00000.Position(3:4) = ps;
                set(fig0 ,'papersize',ps,'paperposition',[0 0 ps])
                set(fig00,'papersize',ps,'paperposition',[0 0 ps])
                set(fig000,'papersize',ps,'paperposition',[0 0 ps])
                set(fig0000,'papersize',ps,'paperposition',[0 0 ps])
                set(fig00000,'papersize',ps,'paperposition',[0 0 ps])
            end
            %
            axisNum = axisNum+1;
            %% plot Umean with depth contours
            figure(fig0)
            eval(sprintf('pos = ppos%d;',axisNum))
            axes('units','centimeters','position',pos);
            imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,Umean), caxis([-0.5 0.5]), colormap(cmocean('balance'))
            hold on,
            contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
            xline((x([iX1 iX2])-50)/(info.xc-50),'--m')
            yline([-2 2],'--g')
            set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',[-5.5 5.5])%0.5*ph/pw*450/info.lc*[-1 1])
            title(title_string,'fontsize',4,'interpreter','latex')
            if axisNum==1
                ylabel('$(y-y_0)/l_c$ [~]','interpreter','latex')
                xlabel('$(x-x_sl)/x_c$ [~]','interpreter','latex')                
            else
                set(gca,'yticklabel',[])
            end
            %
            %% make subplots of each RC transport/max estimate
            figure(fig00)
            axes('units','centimeters','position',pos);
            imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,rms(U,3,'omitnan')), caxis([0 0.5]), colormap(cmocean('amp'))
            hold on,
            contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
            xline((x([iX1 iX2])-50)/(info.xc-50),'--m')
            yline([-2 2],'--g')
            pos = get(gca,'position');
            set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',[-5.5 5.5])%,0.5*ph/pw*450/info.lc*[-1 1])
            title(title_string,'fontsize',4,'interpreter','latex')%            title(RUN,'fontsize',8)
            if axisNum==1
                ylabel('$(y-y_0)/l_c$ [~]','interpreter','latex')
                xlabel('$(x-x_sl)/x_c$ [~]','interpreter','latex')                
            else
                set(gca,'yticklabel',[])
            end
            %
            %% Alongshore velocity
            figure(fig000)
            eval(sprintf('pos = ppos%d;',axisNum))
            axes('units','centimeters','position',pos);
            imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,Vmean), caxis([-0.5 0.5]), colormap(cmocean('balance'))
            hold on,
            contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
            xline( ([x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))]-50)/(info.xc-50),'--c')
            yline([-1 1],'--b')
            yline([-5 5],'--c')
            set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',[-5.5 5.5])%,0.5*ph/pw*450/info.lc*[-1 1])
            %
            title(title_string,'fontsize',4,'interpreter','latex')%            title(RUN,'fontsize',8)
            if axisNum==1
                ylabel('$(y-y_0)/l_c$ [~]','interpreter','latex')
                xlabel('$(x-x_sl)/x_c$ [~]','interpreter','latex')                
            else
                set(gca,'yticklabel',[])
            end
            %
            %% make subplots of mean alongshore current
            figure(fig0000)
            axes('units','centimeters','position',pos);
            imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,rms(V,3,'omitnan')), caxis([0 0.5]), colormap(cmocean('amp'))
            hold on,
            contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5)
            xline( ([x(find(iX0==1,1,'first')) x(find(iX0==1,1,'last'))]-50)/(info.xc-50),'--c')
            yline([-1 1],'--b')
            pos = get(gca,'position');
            set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',[0 450]/(info.xc-50),'ylim',[-5.5 5.5])%,0.5*ph/pw*450/info.lc*[-1 1])
            %
            title(title_string,'fontsize',4,'interpreter','latex')%            title(RUN,'fontsize',8)
            if axisNum==1
                ylabel('$(y-y_0)/l_c$ [~]','interpreter','latex')
                xlabel('$(x-x_sl)/x_c$ [~]','interpreter','latex')                
            else
                set(gca,'yticklabel',[])
            end
            %
            %% Vorticity plot
            figure(fig00000)
            eval(sprintf('pos = ppos%d;',axisNum))
            axes('units','centimeters','position',pos);
            imagesc((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,VORT_mean), caxis([-0.025 0.025]), colormap(cmocean('curl'))
            hold on,
            [~,h] = contour((x-50)/(info.xc-50),(y-info.Ly/2)/info.lc,dep,[0:1:8],'-k','linewidth',0.5,'edgealpha',0.5);
            %% add quiver
            spanx = round(20*info.xc/(200*info.dx));
            spany = round(40*info.lc/(100*info.dy));
            [xx,yy] = meshgrid((x(1:spanx:end)-50)/(info.xc-50),(y(1:spany:end)-info.Ly/2)/info.lc);
            uu = Umean(1:spany:end,1:spanx:end);
            vv = Vmean(1:spany:end,1:spanx:end);
            xlims = [0 450]/(info.xc-50);
            ylims = [-5.5 5.5]; % 0.5*(ph/pw)*(400/info.lc)*[-1 1];
            scalex = 30*(info.xc-50)/info.lc/(info.xc-50);
            scaley = 30*(diff(ylims)/diff(xlims))*(pw/ph)/(info.lc);
            quiver(xx(:),yy(:),uu(:)*scalex, vv(:)*scaley,'ShowArrowHead','off','Marker','.','markersize',1)
            quiver(1.8,4,0.5*scalex, 0,'ShowArrowHead','off','Marker','.','markersize',1,'color','k','linewidth',1)
            text(  1.8,4.75,'0.5 m/s','interpreter','latex','fontsize',5)
            set(gca,'ticklabelinterpreter','latex','fontsize',8,'tickdir','out','ydir','normal','xlim',xlims,'ylim',ylims)
            %
            title(title_string,'fontsize',4,'interpreter','latex')            %            title(RUN,'fontsize',8)
            if axisNum==1
                ylabel('$(y-y_0)/l_c$ [~]','interpreter','latex')
                xlabel('$(x-x_sl)/x_c$ [~]','interpreter','latex')                
            else
                set(gca,'yticklabel',[])
            end
            %
            %% check to see if next grid is the same or if we're at the end of this set of simulations
            if rr==Ndirs
                plotFlag=1;
            else
                nextRUN  = run_dirs{rr+1};
                %% get grid 
                grid_waves = split(nextRUN,'_');
                next_grid  = char(grid_waves{1});
                if ~strcmp(next_grid(1:end-1),grid(1:end-1))
                    plotFlag=1;
                end
            end
            %
            if plotFlag
                figure(fig0);
                cm    = cmocean('balance');
                clims = [-0.5 0.5];
                cvals = clims(1):diff(clims)/255:clims(2);
                cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
                imagesc(0,cvals,reshape(cm,256,1,3))
                set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
                       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
                       'fontsize',6,'tickdir','out')
                xlabel('$\langle u\rangle$ [m/s]','interpreter','latex')
                
                figname = [figDIR,'mean_x_speed_',CASE,'_',RUN,'.png'];
                exportgraphics(fig0,figname)
                
                figure(fig00);
                cm    = cmocean('amp');
                clims = [0 0.5];
                cvals = clims(1):diff(clims)/255:clims(2);
                cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
                imagesc(0,cvals,reshape(cm,256,1,3))
                set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
                       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
                       'fontsize',6,'tickdir','out')
                xlabel('rms$(\bar{u})$ [m/s]','interpreter','latex')
                figname = [figDIR,'rms_x_speed_',CASE,'_',RUN,'.png'];
                exportgraphics(fig00,figname)
                %%
                figure(fig000);
                cm    = cmocean('balance');
                clims = [-0.5 0.5];
                cvals = clims(1):diff(clims)/255:clims(2);
                cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
                imagesc(0,cvals,reshape(cm,256,1,3))
                set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
                       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
                       'fontsize',6,'tickdir','out')
                xlabel('$\langle v\rangle$ [m/s]','interpreter','latex')
                
                figname = [figDIR,'mean_y_speed_',CASE,'_',RUN,'.png'];
                exportgraphics(fig000,figname)
                %
                figure(fig0000);
                cm    = cmocean('amp');
                clims = [0 0.5];
                cvals = clims(1):diff(clims)/255:clims(2);
                cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
                imagesc(0,cvals,reshape(cm,256,1,3))
                set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
                       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
                       'fontsize',6,'tickdir','out')
                xlabel('rms$(\bar{v})$ [m/s]','interpreter','latex')
                figname = [figDIR,'rms_y_speed_',CASE,'_',RUN,'.png'];
                exportgraphics(fig0000,figname)
                %
                
                %
                figure(fig00000);
                cm    = cmocean('curl');
                clims = [-0.025 0.025];
                cvals = clims(1):diff(clims)/255:clims(2);
                cb = axes('units','centimeters','position',cbpos,'ticklabelinterpreter','latex');
                imagesc(0,cvals,reshape(cm,256,1,3))
                set(cb,'ydir','normal','yaxislocation','right','ylim',clims,...
                       'xtick',[],'xaxislocation','top','ticklength',2*get(cb,'ticklength'),...
                       'fontsize',6,'tickdir','out')
                xlabel('$\langle\omega\rangle$ [1/s]','interpreter','latex')
                figname = [figDIR,'vorticity_',CASE,'_',RUN,'.png'];
                exportgraphics(fig00000,figname)
                %
                plotFlag=0;
                axisNum =0;
            end
        end
    end
end
    
    
%% 7) Archive
fileOut = [outDIR,'BulkVelocityStats.mat'];
save(fileOut, 'x',      'y' ,'INFO',...
       'dETA', 'Uscale'     , 'Uscale0',...     
     'Uex','Uex_mean'    ,'Uex_eddy',...
     'Uex_vs_x','Uex_mean_vs_x'    ,'Uex_eddy_vs_x',...     
        'ENS','ENS_mean'    ,'ENS_eddy',...
       'Umax','Umax_mean'   ,'Umax_eddy',...
     'Umax_vs_t','Umax_eddy_vs_t',...
     'Urms','Urms_mean','Urms_eddy',...
     'Umax_BP','Umax_mean_BP','Umax_eddy_BP',...
     'Umax_vs_t_BP','Umax_eddy_vs_t_BP',...     
'Vmean_trough','Veddy_trough','Umean_crest',...
 'Ueddy_crest','Umean_BP'    ,'Ueddy_BP',...
     'Uex_channel', 'Uex_mean_channel' , 'Uex_eddy_channel',...
'Uex_channel_vs_x', 'Uex_mean_channel_vs_x' , 'Uex_eddy_channel_vs_x',...     
     'Uex_ambient', 'Uex_mean_ambient' , 'Uex_eddy_ambient',...
'Uex_ambient_vs_x', 'Uex_mean_ambient_vs_x' , 'Uex_eddy_ambient_vs_x',...     
'Urms_channel', 'Urms_mean_channel' , 'Urms_eddy_channel',...
'Urms_ambient', 'Urms_mean_ambient' , 'Urms_eddy_ambient',...
'ENS_channel', 'ENS_mean_channel' , 'ENS_eddy_channel',...
'ENS_ambient', 'ENS_mean_ambient' , 'ENS_eddy_ambient',...
'height', 'period', 'spread', 'direction',...
'bar_location', 'bar_width', 'bar_amplitude',...
'channel_length', 'channel_amplitude_ratio')
