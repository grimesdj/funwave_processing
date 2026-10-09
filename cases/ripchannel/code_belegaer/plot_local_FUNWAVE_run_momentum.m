clear all
close all
addpath('~/git/funwave/code/')
%
runBATHY = 'spreadRipLong'
runID    = 'barRip0_h10t10s10d00';
%
info = prep_local_ripchannel_info(runBATHY,runID);
%
% info = plot_FUNWAVE_run_momentum(info);

momFile = dir([info.rootMat,'*',runID,'*MomentumTerms.nc']);
momFile = [momFile(1).folder,filesep,momFile(1).name]
info.momFile = momFile;

x = ncread(info.momFile,'x');
y = ncread(info.momFile,'y');

depFile = dir([info.rootMat,'*',runID,'*dep.nc']);
depFile = [depFile(1).folder,filesep,depFile(1).name];
h       = ncread(depFile,'dep');
%        along-shore terms:
%             advection,
FRY  = ncread(momFile,'FRCY');
%
DyVVH = ncread(momFile,'DyVVH');
DxUVH = ncread(momFile,'DxUVH');
ADY  = mean(DyVVH,3,'omitnan') + mean(DxUVH,3,'omitnan'); 
%             pressure grad,
PgrdY = ncread(momFile,'PgrdY');
PGY   = mean(PgrdY,3,'omitnan');
PGYstd= std (PgrdY,[],3,'omitnan');
%             radiation,
DySyy = ncread(momFile,'DySyy');
DxSxy = ncread(momFile,'DxSxy');
BrkDissY = ncread(momFile,'BrkDissY');
RSY     = mean(DySyy,3,'omitnan') + mean(DxSxy,3,'omitnan') - mean(BrkDissY,3,'omitnan');
%
%
%% 2.2) Decompose advection terms into mean and eddy:
ETAmean    = ncread(info.momFile,'etamean');
H      = h+ETAmean;
mask   = (h+ETAmean)>0.1;
%
Umean    = ncread(info.momFile,'umean');
Vmean    = ncread(info.momFile,'vmean');
%
Hmean  = h+ETAmean;
Havg   = mean(Hmean,3,'omitnan');
%
tmp = mean(Umean.*Umean.*Hmean,3,'omitnan');
DxUUHavg = 0*Havg;
DxUUHavg(:,2:end-1) = 0.5*(tmp(:,3:end)-tmp(:,1:end-2))/info.dx;
DxUUHavg(:,[1 end]) = (tmp(:,[2 end])-tmp(:,[1 end-1]))/info.dx;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DyUVHavg = 0*Havg;
DyUVHavg(2:end-1,:) = 0.5*(tmp(3:end,:)-tmp(1:end-2,:))/info.dy;
DyUVHavg([1 end],:) = (tmp([2 end],:)-tmp([1 end-1],:))/info.dy;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DyVVHavg = 0*Havg;
DyVVHavg(2:end-1,:) = 0.5*(tmp(3:end,:)-tmp(1:end-2,:))/info.dy;
DyVVHavg([1 end],:) = (tmp([2 end],:)-tmp([1 end-1],:))/info.dy;
%
tmp = mean(Umean.*Vmean.*Hmean,3,'omitnan');
DxUVHavg = 0*Havg;
DxUVHavg(:,2:end-1) = 0.5*(tmp(:,3:end)-tmp(:,1:end-2))/info.dy;
DxUVHavg(:,[1 end]) = (tmp(:,[2 end])-tmp(:,[1 end-1]))/info.dy;
%
% $$$ ADXavg  = (DxUUHavg+DyUVHavg);
% $$$ ADXeddy = ADX-ADXavg;
ADYavg  = (DyVVHavg+DxUVHavg);
ADYeddy = ADY-ADYavg;
%
%
% only look between shoreline and break-point,
% and +/- 1.5 times the channel width alongshore.
iX = find(x>75 & x<=info.xc);
iY = find(y>info.Ly/2-1.5*info.lc & y<info.Ly/2+1.5*info.lc);
%
PGRS = PGY(:,iX)+RSY(:,iX);
in1 = PGRS  - mean(PGRS,1);
in2 = ADYavg(:,iX)- mean(ADYavg(:,iX),1);
[coh_avg, ky, out1, out2, out12] = alongshore_coherence_estimate(info, in1 , in2);
coh_avg = mean(conj(out12).*out12,2)./( mean(out1,2) .* mean(out2,2) );
%
in2 = ADYeddy(:,iX)- mean(ADYeddy(:,iX),1);
[coh_eddy, ky, out1, out2, out12] = alongshore_coherence_estimate(info, in1 , in2);
coh_eddy = mean(conj(out12).*out12,2)./( mean(out1,2) .* mean(out2,2) );
figure, semilogx(ky, coh_avg,'-b',ky,coh_eddy,'-r')
%
%
% $$$ rcPGRS    = PGY(iY,iX) + RSY(iY,iX);
% $$$ rcADYavg  = ADYavg(iY,iX);
% $$$ rcADYeddy  = ADYeddy(iY,iX);

cm = cmocean('balance');
figure, imagesc(x,y,PGY), colormap(cm), caxis([-1 1]*1e-3)
figure, imagesc(x,y,PGY+RSY), colormap(cm), caxis([-1 1]*1e-3)
figure, imagesc(x,y,PGY+RSY+ADYavg), colormap(cm), caxis([-1 1]*1e-3)
figure, imagesc(x,y,PGY+RSY+ADYavg+ADYeddy), colormap(cm), caxis([-1 1]*1e-3)
figure, imagesc(x,y,PGY+RSY+ADYavg+ADYeddy), colormap(cm), caxis([-1 1]*1e-3)


figure, plot( x, rms(PGY,1),'-k', x, rms(PGY+RSY,1),'-b', x, rms(PGY+RSY+ADY),'-r')
% $$$ figure, imagesc(y,x,PSIavg'), colormap(cm), colorbar, caxis([-0.0118    0.0118]),set(gca,'ydir','normal')


