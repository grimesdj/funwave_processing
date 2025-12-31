function info = estimate_FUNWAVE_run_statistics_WaveAvgVelocity(info);
%
% usage: info = estimate_FUNWAVE_run_statistics_WaveAvgVelocity(info);
%
% calculate run slow time statistics:
% 2) load (u_wavg,v_wavg): 
%    2.1) estimate rotational decomposition (Urot,Vrot), and vorticity
%
% 3) Archive wave average (30-seconds) fields: 
%    (Urot,Vrot, VORT), Fbr, eta, and above statistics and archive 

%
% full alongshore domain
if ~isfield(info,'subDomain')
    iX = find(x0>=0  & x0<=400);
    sdx = length(iX);
    sdy = length(y0);
    subDomain = [1 sdy iX(1) iX(end)];
else
    subDomain = info.subDomain;
end
%
fin = [info.rootMat,info.rootName,'dep.nc'];
h = ncread(fin,'dep',[subDomain([1 3])] , [subDomain([2 4])]);
x = ncread(fin, 'x' , subDomain([3])    ,  subDomain([4]))';
y = ncread(fin, 'y' , subDomain([1])    ,  subDomain([2]));
% $$$ x = x0(subDomain(3):subDomain(4));
% $$$ y = y0(subDomain(1):subDomain(2));
% $$$ h = h0(subDomain(1):subDomain(2), subDomain(3):subDomain(4));
nx = length(x);
ny = length(y);
dx = x(2)-x(1);
dy = y(2)-y(1);
%
% pre-define file names
info.rotVelFile          = [info.rootMat,info.rootName,'velocity_decomposition.nc'];
%
% get sea-surface files:
eta_wavg_files  = dir([info.rootMat,info.rootName,'etawavg*.nc']);
% remove the time averages from the momentum files
momFile  = dir([info.rootMat,info.rootName,'MomentumTerms.nc']);
momFile  = [momFile(1).folder, filesep, momFile.name];
ETAmean = mean(ncread(momFile,'etamean'), 3, 'omitnan');
Umean =   mean(ncread(momFile,'umean'  ), 3, 'omitnan');
Vmean =   mean(ncread(momFile,'vmean'  ), 3, 'omitnan');
%
%
%
%% Calculate rotational current decomposition
files  = dir([info.rootMat,info.rootName,'uwavg*.nc']);
Nf     = length(files);
t0     = [];
Urot   = [];
Vrot   = [];
PSI    = [];
VORT   = [];
ERR    = [];
eta0   = [];
iter = 1;
for ii=1:Nf
    %
    fprintf('loading from: %s \n', files(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'uwavg_%02d.nc'],ii);
    t   = ncread(fin,'t');
    nt  = length(t);
    %
    % 3) vorticity statistics
    vars = {'etawavg','uwavg','vwavg'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
    end
    %
    eta=etawavg-ETAmean;
    u = uwavg-Umean;
    v = vwavg-Vmean;
    clear etawavg uwavg vwavg
    %
    % estimate vorticity
    [uy,~ ,~] = gradientDG(u./dy);
    [~ ,vx,~] = gradientDG(v./dx);
    omega = vx-uy;
    clear uy vx
    %
    %
    %
    % helmholtz decomposition
    II=0;
    for jj = 1:nt
        II=II+1;
        U = u(:,:,jj);
        V = v(:,:,jj);
        [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(U,V,dx,dy);
        urot(:,:,II) = u_psi;
        vrot(:,:,II) = v_psi;
        psi0(:,:,II) = psi;
        err (:,:,II) = sqrt( (U-(u_psi+u_phi)).^2 + (V-(v_psi+v_phi)).^2 );
    end
    % log Urot, Vrot
    t0   = cat(1,t0,t);
    VORT = cat(3,VORT,omega);
    Urot = cat(3,Urot,urot);
    Vrot = cat(3,Vrot,vrot);
    PSI = cat(3,PSI,psi0);
    ERR = cat(3,ERR,err);
    eta0= cat(3,eta0,eta);
    clear err psi0 urot vrot u_psi v_psi U V 
    %
    %
end
%
save(info.fileName,'-struct','info')
%
% For debugging the netcdf write portion, uncomment this:
% $$$ save('/scratch/grimesdj/ripchannel/planar2D/mat_data/debugging_code.mat','-v7.3')
%
% 2) info.rotVelFile          = [info.rootMat,info.rootName,'velocity_decomposition.nc'];
% 2.1) variables:
%      VORT = low-pass vorticity
%      Urot = low-pass rotational cross-shore velocity
%      Vrot = low-pass rotational alongshore velocity
%      PSI  = low-pass streamfunction
%      ERR  = error between decomposition and actual
if exist(info.rotVelFile,'file')
    eval(['!rm ',info.rotVelFile])
end
% this is for archiving
dim_yx  = {"y",length(y),"x",length(x)};
dim_yxt = {"y",length(y),"x",length(x),"t",length(t0)};
%
nccreate  (info.rotVelFile,'VORT','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'VORT',VORT);
ncwriteatt(info.rotVelFile,'VORT','Description','low-pass vorticity');
%
nccreate  (info.rotVelFile,'x','Dimensions',{"x",length(x)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'x',x);
nccreate  (info.rotVelFile,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'y',y);
nccreate  (info.rotVelFile,'t','Dimensions',{"t",length(t0)},'Format','netcdf4')
ncwrite   (info.rotVelFile,'t',t0);
%
nccreate  (info.rotVelFile,'Urot','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'Urot',Urot);
ncwriteatt(info.rotVelFile,'Urot','Description','low-pass rotational cross-shore velocity');
%
nccreate  (info.rotVelFile,'Vrot','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'Vrot',Vrot);
ncwriteatt(info.rotVelFile,'Vrot','Description','low-pass rotational alongshore velocity');
%
nccreate  (info.rotVelFile,'PSI','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'PSI',PSI);
ncwriteatt(info.rotVelFile,'PSI','Description','low-pass streamfunction');
%
nccreate  (info.rotVelFile,'ERR','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'ERR',ERR);
ncwriteatt(info.rotVelFile,'ERR','Description','low-pass error in helmholtz decomposition');
%
nccreate  (info.rotVelFile,'eta','Dimensions',dim_yxt,'Format','netcdf4')
ncwrite   (info.rotVelFile,'eta',eta0);
ncwriteatt(info.rotVelFile,'eta','Description','low-pass sea-surface elevation');
%
%
return
% $$$ 
% $$$ 
% $$$ 
% $$$ % this if from the wave stats code
% $$$ eta = eta0;
% $$$ t   = t0;
% $$$ nt  = length(t);
% $$$ dt  = mean(diff(t));
% $$$ %
% $$$ [Hs_xy,eta_bar,mask0,x,y,h,freq,Snn_xy,xsl] = calculate_funwave_wave_height_statistics_v2(info.rootMat,info.rootName,info.bathyFile,info.Hs,info.Tp,subDomain);
% $$$ %
% $$$ info.x_shoreline = xsl;
% $$$ info.mask        = mask0;
% $$$ info.eta_bar     = eta_bar;
% $$$ %
% $$$ fig0 = figure;
% $$$ Hs_x = nanmean(Hs_xy,1);
% $$$ p0 = plot(x,Hs_xy,'.k',x,Hs_x,'-r','markersize',1,'linewidth',1.5);
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$H_s$ [m]','interpreter','latex')
% $$$ set(gca,'xlim',[75 600],'ylim',[0 1.1*nanmax(Hs_xy(:))],'ticklabelinterpreter','latex','tickdir','out')
% $$$ f0l1 = legend([p0(1),p0(end)]','$H_\mathrm{s}(x,y)$','$\bar{H}_\mathrm{s}(x)$');
% $$$ set(f0l1,'location','southeast','interpreter','latex')
% $$$ title(info.runName)
% $$$ if ~exist([info.rootSim,filesep,'figures'],'dir')
% $$$     eval(['!mkdir ',[info.rootSim,filesep,'figures']])
% $$$ end
% $$$ exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'Hs.pdf'])
% $$$ %
% $$$ end
% $$$ %
% $$$ % get front statistics
% $$$ [mean_stats,binned_stats] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ % $$$ [mL,sL,N,mLx,sLx,Lx_log_mean,Lx_log_std,Nx,xylog,Xbins] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ % front_file = [info.rootMat,filesep,info.rootName,'bore_front_statistics.mat'];
% $$$ %
% $$$ % make a stats plot? what stats?
% $$$ fig1 = figure;
% $$$ p1 = plot(binned_stats.Xbins,binned_stats.Length,'-k',binned_stats.Xbins,exp( binned_stats.log_mean_length' + binned_stats.log_std_length'*[-1 1] ), '--r');
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$L$ [m]','interpreter','latex')
% $$$ set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
% $$$ f1l1 = legend([p1(1), p1(2)]','$\bar{L}(l)=\exp\left(\overline{\log(l)}\right)$','$\bar{L}\pm\mathrm{std}(L)$');
% $$$ set(f1l1,'location','northeast','interpreter','latex')
% $$$ title(info.runName)
% $$$ if ~exist([info.rootSim,filesep,'figures'],'dir')
% $$$     eval(['!mkdir ',[info.rootSim,filesep,'figures']])
% $$$ end
% $$$ exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'crest_length_vs_x.pdf'])
% $$$ %
% $$$ fig2 = figure;
% $$$ p2 = plot(binned_stats.Xbins,binned_stats.N,'-k','linewidth',2);
% $$$ xlabel('crosshore [m]','interpreter','latex')
% $$$ ylabel('$N$ [crests/frame]','interpreter','latex')
% $$$ set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
% $$$ title(info.runName)
% $$$ exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'crests_per_frame_vs_x.pdf'])
% $$$ %
% $$$ close all
% $$$ %
% $$$ wave_file = [info.rootMat,filesep,info.rootName,'wave_statistics.mat'];
% $$$ info.waveStatsFile = wave_file;
% $$$ if exist(wave_file,'file')
% $$$     save(wave_file,'-append')
% $$$ else
% $$$      save(wave_file,'-v7.3')
% $$$ end
% $$$ save(info.fileName,'-struct','info')
% $$$ %
