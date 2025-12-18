% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% code to be launched on cms-hpc "cuttlefish"
% 0) requires the input bathymetry name as top-dir
runBATHYs = {'test2D','testRip'};
iter = 1;
run_log   = {};
for ii=1:length(runBATHYs)
    runBATHY = runBATHYs{ii};
    %
    runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
    matDIR   = [runDIR,filesep,'mat_data'];
    figDIR   = [runDIR,filesep,'figures/'];
    %
    % the list of run directories are saved in:
    load([matDIR,filesep,'runs_to_process.mat'])
    %
    % loop over run_dirs
    Ndirs  = length(run_dirs);
    run_log= cat(1,run_log,run_dirs);
    for jj = 1:Ndirs
        % 1) get current run subdirectory to process:
        runID    = run_dirs{jj};
        fprintf('\n processing: %s %s \n', runBATHY,runID)    
        % 2) get the archived info structure:
        infoFile = dir([matDIR,filesep,'*','info','*',runID,'.mat']);
        info     = load([infoFile(1).folder,filesep,infoFile(1).name]);
        %
        % files needed for subsequent analysis
        velFile  = info.rotVelFile;
        momFile  = [info.rootMat,info.rootName,'MomentumTerms.nc'];
        %
        x   = ncread(velFile,'x');
        y   = ncread(velFile,'y');
        t   = ncread(velFile,'t');
        %
        if ii==1 & jj==1;
            Uex = nan(length(x),Ndirs);
        end
        %
        % load breaking to estimate surfzone width
        BrkDissX = ncread(momFile,'BrkDissX');
        tmp0     = BrkDissX; tmp0(~info.mask)=nan;  rms0 = rms(tmp0,[1 3],'omitnan');
        cum_rms  = cumsum(rms0);
        idx_Lsz  = find(cum_rms>=0.98*max(cum_rms),1,'first');
        Lsz(iter)  = x(idx_Lsz);
        info.Lsz = x(idx_Lsz);
        save(info.fileName,'-struct','info')
        %
        % now get rotational cross-shore velocity to estimate Uex and Hovmoller
        Urot     = ncread(velFile,'Urot');
        HOV      = reshape(Urot(:,idx_Lsz,:),[length(y),length(t)]);
        Uex(:,iter)= sum(Urot.*(Urot>0),[1 3])./sum((Urot>0),[1 3]);
        %
        % estimate the alongshore wavenumber spectra of HOV
        [Suu,ky] = alongshore_spectra_estimate(info,HOV);
        Suu = mean(Suu,2);
        Sex(:,iter)= Suu;
        %
        % make a hovmoller plot
        xm = 2.5;
        ym = 2.5;
        pw = 8;
        ph = 3;
        ag = 0.5;
        ppos1 = [xm ym        pw ph];
        ppos2 = [xm ym+ph+ym  pw ph];
        cbpos = [xm+pw+ag ym ag ph*2/3];
        ps    = [2*xm+pw+3*ag, 3.5*ym+2*ph+ag];
        fig   = figure('units','centimeters');
        fig.Position(3:4)=ps;
        fig.PaperSize=ps;
        fig.PaperPosition=[0 0 ps];
        %
        clims = 1.5*std(HOV(:))*[-1 1];
        cm    = cmocean('balance');
        clrs  = clims(1):diff(clims)/255:clims(2);
        %
        a1 = axes('units','centimeters','position',ppos1);
        imagesc(y,(t-t(1)),HOV')
        colormap(cm), caxis(clims)
        xlabel('$y$ [m]','interpreter','latex')
        ylabel('$t$ [s]','interpreter','latex')
        set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
        %
        a2 = axes('units','centimeters','position',ppos2);
        semilogx(ky,Suu,'-k','linewidth',2)
        xlabel('$k_y$ [m$^-1$]','interpreter','latex')
        ylabel('$S_{u_\psi u_\psi}$ (m/s)$^2$','interpreter','latex')
        set(a2,'tickdir','out','ticklabelinterpreter','latex')
        %
        cb = axes('units','centimeters','position',cbpos);
        imagesc(0,clrs,reshape(cm,256,1,3))
        xlabel('[m/s]','interpreter','latex')
        set(cb,'xtick',[],'xaxislocation','top','yaxislocation','right','ydir','normal',...
               'ticklabelinterpreter','latex')
        %
        figname = [figDIR,info.runName,'_Urot_Hovmoller_and_Spectra.pdf'];
        exportgraphics(fig,figname)
        close(fig)
        %
        iter = iter+1;
    end
end
%
% Uex and spectra plot
xm = 3;
ym = 2.5;
pw = 8;
ph = 3;
ag = 0.5;
ppos1 = [xm ym        pw ph];
ppos2 = [xm ym+ph+ym  pw ph];
ps    = [2*xm+pw, 3.5*ym+2*ph+ag];
fig   = figure('units','centimeters');
fig.Position(3:4)=ps;
fig.PaperSize=ps;
fig.PaperPosition=[0 0 ps];
%
clims = 1.5*std(HOV(:))*[-1 1];
cm    = cmocean('balance');
clrs  = clims(1):diff(clims)/255:clims(2);
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x./Lsz,100*Uex,'linewidth',2)
colormap(cm), caxis(clims)
ylabel('$U_\mathrm{ex}$ [cm/s]','interpreter','latex')
xlabel('$x/L_\mathrm{sz}$ [s]','interpreter','latex')
set(a1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal')
%
a2 = axes('units','centimeters','position',ppos2);
semilogx(ky,Sex,'-','linewidth',2)
xlabel('$k_y$ [m$^-1$]','interpreter','latex')
ylabel('$S_{u_\psi u_\psi}$ (m/s)$^2$','interpreter','latex')
runNames = regexp(run_log,'^[^_]+(?=_)','match');
runNames = vertcat(runNames{:})';
legend(runNames,'interpreter','latex');
set(a2,'tickdir','out','ticklabelinterpreter','latex')
%
figname = [figDIR,'Uex_and_Spectra.pdf'];
exportgraphics(fig,figname)
close(fig)
