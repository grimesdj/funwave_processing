% code to be launched on cms-hpc "cuttlefish"
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
% requires the input bathymetry name as top-dir
runBATHY = 'resolution1D';%'planar1D'
figDIR = '/storage/cms/grimesdj_lab/grimesdj/git/funwave/cases/ripchannel/testCases/test_1D_resolution/figures/';
%
% reprocess raw output?
reproc=0;
% remove raw text output?
rmfiles=0;
% number of samples in each data-file
Ns = 5500;
%
%
T_INTV_mean = 100;
% number of time-steps to average for breaking wave dissipation rate,
% should be the same as T_INTV_mean for radiation stress calculations:
Navg = T_INTV_mean;
runDIR   = ['/scratch/grimesdj/ripchannel/',runBATHY];
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
%
% loop over run_dirs
Ndirs  = length(run_dirs);
%
% create 2-panel plot for momentum terms
xm = 2.5;
ym = 2;
pw = 8;
ph = 4;
ppos1 = [xm ym        pw ph];
ppos2 = [xm ym+ph+0.5 pw ph];
ps    = [2*xm+pw 1.5*ym+2*ph];
%
fig0   = figure('units','centimeters');
fig0.Position(3:4) = ps;
fig0.PaperSize = ps;
fig0.PaperPosition = [0 0 ps];
%
%
res1D = struct([]);
for jj=1:Ndirs
    % 0) load info file
    runNAME  = run_dirs{jj};
    rawDIR   = [runDIR,filesep,runNAME,filesep,'output'];
    fprintf('\nworking on: %s\n', runNAME)
    infoFile = [matDIR,filesep,'ripchannel_run_info_',runNAME,'.mat'];
    [info,t0,dt0] = ripchannel_run_info_cmshpc(matDIR,runNAME);
    % 1) convert output to netcdf
    if reproc
    spanx = 1;
    spany = 1;
    %
    rngx  = [1 400/info.dx];
    rngy  = [1 info.Ly/info.dy];
    %        
    rootOut = split(info.rootOut,':');
    vars = {'dep','eta','u','v','p','q','mask','BrkSrcX'};
    fLog = convert_funwave_output_to_NetCDF(rootOut{end},[info.rootMat,info.rootName],vars,t0,dt0,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,Ns);
    %
    BrkDissFiles = dir([rootOut{end},filesep,'BrkDissX*']);
    Nbrk = length(BrkDissFiles);
    dt_lp = T_INTV_mean*ones(Nbrk,1);
    t_lp = [1:Nbrk]*dt_lp(1);
    vars = {'BrkDissX','DxSxx','DxUUH','FRCX','PgrdX','Sxx'};
    fLog = convert_funwave_output_to_single_NetCDF(rootOut{end},[info.rootMat,info.rootName,'MomentumTerms'],vars,t_lp,dt_lp,info.dx,spanx,rngx,info.dy,spany,rngy,rmfiles,round(Ns/Navg));
    end
    %
    % need depth, ideally (H+eta_bar)
    depFiles = dir([info.rootMat,info.rootName,'dep*']);
    depFile  = [depFiles(1).folder,filesep,depFiles(1).name];
    H = ncread(depFile,'dep');
    x = ncread(depFile,'x');
    y = ncread(depFile,'y');
    %
    % load and time-average the sea-surface elevation...
    etaFiles = dir([info.rootMat,info.rootName,'eta*']);
    etaFile  = [etaFiles(1).folder,filesep,etaFiles(1).name];
    eta      = ncread(etaFile,'eta');
    H        = H + mean(eta,3);
    Hs       = 4*std(eta,[],3);
    %
    % 2) time-average BrkSrcX to match BrkDissX
    brksrcFiles = dir([info.rootMat,info.rootName,'BrkSrcX*']);
    brksrcFile  = [brksrcFiles(1).folder,filesep,brksrcFiles(1).name];
    BrkSrcX     = ncread(brksrcFile,'BrkSrcX');
    ts          = ncread(brksrcFile,'t');
    ts          = ts-ts(1)+median(diff(ts));
    %
    [Ny,Nx,Nsrc]  = size(BrkSrcX);
    % number of averages available:
    Na = floor(Nsrc/Navg);
    tmp    = reshape(BrkSrcX(:,:,1:Navg*Na),[Ny,Nx,Navg,Na]);
    tmp    = sum(tmp,3)/Navg;
    BrkSrcX = reshape(tmp,[Ny,Nx,Na]);
    %
    tmp    = reshape(ts(1:Navg*Na),[Navg,Na]);
    tmp    = round(tmp(end,:));
    %
    % 3) load momentum terms DxSxx, dep, Sxx
    mtermFiles = dir([info.rootMat,info.rootName,'MomentumTerms*']);
    mtermFile  = [mtermFiles(1).folder,filesep,mtermFiles(1).name];
    tm         = ncread(mtermFile,'t');
    BrkDissX   = ncread(mtermFile,'BrkDissX');    
    DxSxx      = ncread(mtermFile,'DxSxx');
    Sxx        = ncread(mtermFile,'Sxx');    
    DxUUH      = ncread(mtermFile,'DxUUH');
    PgrdX      = ncread(mtermFile,'PgrdX');
    FRCX       = ncread(mtermFile,'FRCX');
    %
    clf(fig0)
    % 4) plot the time-averaged momentum terms (panel 1)
    BrkSrcX_avg  = mean(BrkSrcX,3);
    BrkDissX_avg = mean(BrkDissX,3);
    DxUUH_avg    = mean(DxUUH,3);    
    DxSxx_avg    = mean(DxSxx,3);
    Sxx_avg      = mean(Sxx,3);    
    PgrdX_avg    = mean(PgrdX,3);
    FRCX_avg     = mean(FRCX,3);
    %
    % apply a 5-m filter to all fields to dealias 1Hz output (this should not be necessary for momentum terms!)
    dx = x(2)-x(1);
    Nf = round(20/dx); if ~mod(Nf,2), Nf=Nf+1; end
    flt= hamming(Nf); flt = flt./sum(flt);
    %
    Hs            = conv(Hs,flt,'same');    
    BrkSrcX_avg   = conv(BrkSrcX_avg,flt,'same');
    BrkDissX_avg  = conv(BrkDissX_avg,flt,'same');
    DxUUH_avg     = conv(DxUUH_avg,flt,'same');    
    DxSxx_avg     = conv(DxSxx_avg,flt,'same');
    Sxx_avg       = conv(Sxx_avg,flt,'same');    
    PgrdX_avg     = conv(PgrdX_avg,flt,'same');
    FRCX_avg      = conv(FRCX_avg,flt,'same');
    %
    % Nan values where mean water depth is below 0.1m
    dry               = H<0.1 | x(:)'>x(end-Nf);
    Hs(dry)           = nan;
    BrkSrcX_avg(dry)  = nan;
    BrkDissX_avg(dry) = nan;
    DxSxx_avg(dry)    = nan;
    Sxx_avg(dry)      = nan;    
    DxUUH_avg(dry)    = nan;    
    PgrdX_avg(dry)    = nan;
    FRCX_avg(dry)     = nan;
    %
    Rx  = (PgrdX_avg + DxSxx_avg + FRCX_avg + DxUUH_avg-BrkDissX_avg);
    ax1 = axes('units','centimeters','position',ppos2);
    plot(x,PgrdX_avg,'-k',x,DxSxx_avg,'-r',x,FRCX_avg,'-b',x,DxUUH_avg,'-g',x,-BrkDissX_avg,'-m',x,Rx,'--c','linewidth',1)
    grid on
    set(ax1,'xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    ylabel(ax1,'[m$^2$/s$^2$]','interpreter','latex')
    lg1 = legend('$\partial\eta/\partial x$','$\partial S_{xx}/\partial x$','$\tau_b$','$\partial (u^2 H)/\partial x$','$F_{\mathrm{br},x}$');
    set(lg1,'interpreter','latex','fontsize',8,'location','southeast','color','none')
    %
% $$$     % 5) plot the FRC, SRC, DISS (panel 2)
% $$$     DxSxx_on_H = DxSxx_avg./max(H,0.1);
    % approximate d/dx (Sxx) ~ 3/2 d/dx (E) ~ 3/2 d/dx (Hs^2)
    DxSxx_from_Hs = -3/2*Hs.*gradient(Hs)/dx;
    ax2 = axes('units','centimeters','position',ppos1);
    plot(x,-DxSxx_avg,'-k',x,BrkDissX_avg,'-r',x,BrkSrcX_avg,'-b',x,DxSxx_from_Hs,'-g','linewidth',1)
    grid on
    set(ax2,'tickdir','out','ticklabelinterpreter','latex','fontsize',15)
    set(ax2.YAxis,'Exponent',0)
    ylabel(ax2,'[m$^2$/s$^2$]','interpreter','latex')
    xlabel(ax2,'$x$ [m]','interpreter','latex')    
    lg2 = legend('$-\partial S_{xx}/\partial x$','$\bar{F}_\mathrm{br}$','$\langle F_\mathrm{br}\rangle$','$3/2\partial(H_\mathrm{rms}^2)/\partial x$');
    set(lg2,'interpreter','latex','fontsize',8,'location','southeast','color','none')
    figname = [figDIR,'momentum_terms_',runNAME,'.pdf'];
    exportgraphics(fig0,figname)
    %
    % 6) log time-averaged quantities for archiving
    res1D(jj).x = x;
    res1D(jj).H = H;
    res1D(jj).Hs = Hs;    
    res1D(jj).BrkSrcX = BrkSrcX_avg;
    res1D(jj).BrkDissX= BrkDissX_avg;
    res1D(jj).DxSxx   = DxSxx_avg;
    res1D(jj).Sxx     = Sxx_avg;    
    res1D(jj).PgrdX   = PgrdX_avg;
    res1D(jj).FRCX    = FRCX_avg;
    res1D(jj).DxUUH   = DxUUH_avg;
    %
end
%
% 7) archive time-averages
archive = [matDIR,filesep,runBATHY,'_time_averaged_momentum_terms.mat'];
save(archive,'res1D')
