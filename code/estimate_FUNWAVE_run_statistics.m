function info = estimate_FUNWAVE_run_statistics(info);
%
% usage: info = estimate_FUNWAVE_run_statistics(info);
%
% calculate run fast/slow time statistics:
% 1) load (eta,mask,nubrk): save these in structure... keep until end!
%    1.1) wave frequency spectra, wave height,
%    1.2) breking front statistics
% 2) load (u,v,p,q): average these as we go... save a padded region from end of each
%    2.1) estimate rotational decomposition (Urot,Vrot), and vorticity
%    2.2) estimate breaking dissipation rate and wave force
%    2.3) estimate budget terms:
%    2.3.1) Rotational power: <u_rot dot Fbr>, <u_rot> dot <Fbr>, <u_rot' dot Fbr'>
%                             coh(u_rot,Fbr), 
%    2.3.1) Rotational impulse: <vort dot Fbr>, <vort> dot <Fbr>, <vort' dot Fbr'>
%                              coh(vort,Fbr), coh( int(vort), int(Fbr) )
%
% 3) wave average (30-seconds) fields: 
%    (Urot,Vrot, VORT), Fbr, eta, and above statistics and archive 

%
% load the depth and ancillary fields
depFile = [info.rootMat,info.rootName,'dep.nc'];
t0 = ncread(depFile,'t');
nt0= length(t0);
x0 = ncread(depFile,'x');
y0 = ncread(depFile,'y');
h0 = ncread(depFile,'dep');
dt = t0(2)-t0(1);
dx = x0(2)-x0(1);
dy = y0(2)-y0(1);
[xx,yy] = meshgrid(x0,y0);
%
% define the filter parameters:
% filter halfwidth:
Ns  = round(30/dt); if mod(nf,2), nf=nf+1; end
% filter width
nf  = 2*Ns+1;
flt = hamming(nf); flt=flt/sum(flt);
%
% full alongshore domain
if ~isfield(info,'subDomain')
    iX = find(x0>=0  & x0<=400);
    nx = length(iX);
    ny = length(y0);
    subDomain = [1 ny iX(1) iX(end)];
else
    subDomain = info.subDomain;
end
%
x = x0(subDomain(2):subDomain(4));
y = y0(subDomain(1):subDomain(3));
h = h0(subDomain(1):subDomain(2), subDomain(3):subDomain(4));
%
% pre-define file names
info.waveForceFile       = [info.rootMat,info.rootName,'wave_forcing.nc'];
info.rotVelFile          = [info.rootMat,info.rootName,'velocity_decomposition.nc'];
info.waveStatsFile       = [info.rootMat,info.rootName,'wave_statistics.nc'];
%
%
% load data
files  = dir([info.rootMat,info.rootName,'eta*.nc']);
Nf     = length(files);
nubrk0 = [];
eta0   = [];
t0     = [];
%
% preallocate pads:
p_pad = []; q_pad = []; Fx_pad=[]; Fy_pad=[];
% preallocate archive vars:
Eb0=[]; Ebavg=[]; Ib0 = []; VORT = []; Urot = []; Vrot = [];PSI=[];
for ii=1:Nf
    fprintf('loading eta from: %s \n', files(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'eta_%02d.nc'],ii);
    ncid= netcdf.open(fin);
    [~,nt]  = netcdf.inqDim(ncid,2);
    netcdf.close(ncid)
    eta     = ncread(fin,'eta',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);
    t       = ncread(fin,'t');
    %
    % 
    % 2) estimate wave forcing, power, vorticity
    % 2.0) read in (p,q,nubrk), estimate Fbr,
    vars = {'p','q','nubrk'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
    end
    H = h+eta;
    % 2.2) estimate breaking dissipation rate and wave force
    [Fbx, Fby] = estimate_Fbr(p,q,nubrk,H,dx,dy);
    %
    % 2) wave average forcing statistics
    p = cat(3,p_pad, p);
    q = cat(3,q_pad, q);
    Fbx = cat(3,Fx_pad, Fbx);
    Fby = cat(3,Fy_pad, Fby);
    Nt  = size(Fby,3);
    % 
    % fractional number of full filter segements
    ns = Nt/Ns;
    % number of whole segments
    N  = floor( ns );
    % remainder
    ns = rem(ns,1);   
    %
    pavg = time_average_field(p, flt);
    qavg = time_average_field(q, flt);    
    Eb   = Fbx .* (p-pavg) + Fby .* (q-qavg);
    Eb   = time_average_field(Eb, flt);
    %
    Eb0   = cat(3,Eb0   ,Eb   (:,:,(Ns:Ns:Ns*(N-1))));
    eta0  = cat(3,eta0  ,eta  (:,:,(Ns:Ns:Ns*(N-1))));
    nubrk0= cat(3,nubrk0,nubrk(:,:,(Ns:Ns:Ns*(N-1))));              
    % keep end points
    np    = size(p_pad,3);
    idt   = round((np-1)+Nt-Ns*(1+ns)+1:Nt); 
    p_pad = p(:,:,idt);
    q_pad = q(:,:,idt);
    clear p q Eb
    %
    % estimate mean product
    Fbx_avg = time_average_field(Fbx, flt);
    Fby_avg = time_average_field(Fby, flt);
    %
    tmp      = Fbx_avg .* p_avg + Fby_avg .* q_avg;
    Ebavg    = cat(3,Ebavg   ,tmp(:,:,(Ns:Ns:Ns*(N-1))));
    %
    % 3) vorticity statistics
    vars = {'u','v'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''',[subDomain([1 3]), 1] , [subDomain([2 4]), nt]);'])
    end
    %
    p = cat(3,p_pad, p);
    q = cat(3,q_pad, q);
    %
    % estimate vorticity
    [uy,~ ,~] = gradientDG(u./dy);
    [~ ,vx,~] = gradientDG(v./dx);
    omega = vx-uy;
    clear uy vx
    %
    [~,~,dO] = gradientDG(omega./dt);
    %
    % estimate curl of Fbr
    [Fbx_y,  ~  ] = gradientDG(Fbx./dy);
    [~    ,Fby_x] = gradientDG(Fby./dx);        
    cFbr  = Fby_x - Fbx_y;
    clear Fbx_y Fby_x
    cI    = cumsum(cFbr*dt,3);
    %
    % estimate statistics
    tmp = struct('dy',dy,'Ny',length(y0));
    [~,ky,dO_spec,cI_spec,cO_cI_cospec] = alongshore_coherence_estimate(tmp,dO,cI);
    clear dO cI
    %
    % wave average spectra
    coh = (cO_cI_cospec)./sqrt( cO_spec.*cI_spec );
    phi = atan2(imag(coh),real(coh));
    coh = abs(coh).^2;
    coh = time_average_field(coh, flt);
    %
    Ib0   = cat(3,Ib0   ,coh   (:,:,(Ns:Ns:Ns*(N-1))));
    clear coh phi dO_spec cI_spec cO_cI_cospec
    %
    %
    Fx_pad= Fbx(:,:,idt);
    Fy_pad= Fby(:,:,idt);
    clear Fbx Fby
    %
    % helmholtz decomposition 
    % first time average
    u_avg = time_average_field(u, flt);
    v_avg = time_average_field(v, flt);
    u_avg = u_avg(:,:,(Ns:Ns:Ns*(N-1)));
    v_avg = v_avg(:,:,(Ns:Ns:Ns*(N-1)));
    II=0;
    for jj = 1:N*Ns
        II=II+1;
        U = u_avg(:,:,jj);
        V = v_avg(:,:,jj);
        [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(U,V,dx,dy);
        urot(:,:,II) = u_psi;
        vrot(:,:,II) = v_psi;
        psi0(:,:,II) = psi;
        err (:,:,II) = sqrt( (U-(u_psi+u_phi)).^2 + (V-(v_psi+v_phi)).^2 );
    end
    % log Urot, Vrot
    VORT = cat(3,VORT,omega(:,:,(Ns:Ns:Ns*(N-1))));
    Urot = cat(3,Urot,urot);
    Vrot = cat(3,Vrot,vrot);
    PSI = cat(3,PSI,psi0);
    ERR = cat(3,ERR,err);                
    clear err psi0 urot vrot u_psi v_psi U V u_avg v_avg
    %
    u_pad = u(:,:,idt);
    v_pad = v(:,:,idt);
    clear u v
    %
    % 4) archive eta and nubrk
    eta0 = cat(3,eta0,eta);
    nubrk0= cat(3,nubrk0,nubrk);
    t0   = cat(1,t0,t);
    clear eta nubrk t
end
save('/scratch/grimesdj/ripchannel/planar2D/mat_data/debugging_code.mat')
return
% this is for archiving
nccreate(fout,var,'Dimensions',dim,'Format','netcdf4')
eval(['ncwrite(fout,''',var,''',',var,');'])
% log this averaged field
dim   = {"y",length(y),"x",length(x),"t",length(t)};


% this if from the wave stats code
eta = eta0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
%
[Hs_xy,eta_bar,mask0,x,y,h,freq,Snn_xy,xsl] = calculate_funwave_wave_height_statistics_v2(info.rootMat,info.rootName,info.bathyFile,info.Hs,info.Tp,subDomain);
%
info.x_shoreline = xsl;
info.mask        = mask0;
info.eta_bar     = eta_bar;
%
fig0 = figure;
Hs_x = nanmean(Hs_xy,1);
p0 = plot(x,Hs_xy,'.k',x,Hs_x,'-r','markersize',1,'linewidth',1.5);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$H_s$ [m]','interpreter','latex')
set(gca,'xlim',[75 600],'ylim',[0 1.1*nanmax(Hs_xy(:))],'ticklabelinterpreter','latex','tickdir','out')
f0l1 = legend([p0(1),p0(end)]','$H_\mathrm{s}(x,y)$','$\bar{H}_\mathrm{s}(x)$');
set(f0l1,'location','southeast','interpreter','latex')
title(info.runName)
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ',[info.rootSim,filesep,'figures']])
end
exportgraphics(fig0,[info.rootSim,filesep,'figures',filesep,info.rootName,'Hs.pdf'])
%
end
%
% get front statistics
[mean_stats,binned_stats] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% $$$ [mL,sL,N,mLx,sLx,Lx_log_mean,Lx_log_std,Nx,xylog,Xbins] = compile_funwave_bore_fronts(info.rootMat,info.rootName,info.bathyFile,subDomain);
% front_file = [info.rootMat,filesep,info.rootName,'bore_front_statistics.mat'];
%
% make a stats plot? what stats?
fig1 = figure;
p1 = plot(binned_stats.Xbins,binned_stats.Length,'-k',binned_stats.Xbins,exp( binned_stats.log_mean_length' + binned_stats.log_std_length'*[-1 1] ), '--r');
xlabel('crosshore [m]','interpreter','latex')
ylabel('$L$ [m]','interpreter','latex')
set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
f1l1 = legend([p1(1), p1(2)]','$\bar{L}(l)=\exp\left(\overline{\log(l)}\right)$','$\bar{L}\pm\mathrm{std}(L)$');
set(f1l1,'location','northeast','interpreter','latex')
title(info.runName)
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ',[info.rootSim,filesep,'figures']])
end
exportgraphics(fig1,[info.rootSim,filesep,'figures',filesep,info.rootName,'crest_length_vs_x.pdf'])
%
fig2 = figure;
p2 = plot(binned_stats.Xbins,binned_stats.N,'-k','linewidth',2);
xlabel('crosshore [m]','interpreter','latex')
ylabel('$N$ [crests/frame]','interpreter','latex')
set(gca,'xlim',[75 300],'ticklabelinterpreter','latex','tickdir','out')
title(info.runName)
exportgraphics(fig2,[info.rootSim,filesep,'figures',filesep,info.rootName,'crests_per_frame_vs_x.pdf'])
%
close all
%
wave_file = [info.rootMat,filesep,info.rootName,'wave_statistics.mat'];
info.waveStatsFile = wave_file;
if exist(wave_file,'file')
    save(wave_file,'-append')
else
     save(wave_file,'-v7.3')
end
save(info.fileName,'-struct','info')
%
