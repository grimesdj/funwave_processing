function fout = plot_FUNWAVE_offline_momentum_budget(info);
%
% USAGE: fout = plot_FUNWAVE_offline_momentum_budget(info);
%
% Use offline wave averaged velocities (U+u), 
% where \hat{u} = U + u + \tilde{u}, with \tilde{u}=wave-velocity,
% to estimate momentum budget terms:
% 
% crossshore (x): d/dx( UU + uu ) + d/dy (UV + uv)
% alongshore (y): d/dx( UV + uv ) + d/dy (VV + vv)
%
% compare these to the online calculated pressure gradients:
%
% cross-shore (x): -g d/dx \eta
% along-shore (y): -g d/dy \eta

% momentum file
momFile = [info.rootMat,info.rootName,'MomentumTerms.nc'];
%
% $$$ %
% load from the momentum file:
x   = ncread(momFile,'x');
y   = ncread(momFile,'y');
t   = ncread(momFile,'t');
%
% load the depth for 2D plots
depFile = [info.rootMat,info.rootName,'dep.nc'];
h = ncread(depFile,'dep');
%
%
ETA   = ncread(momFile,'etamean');
U     = ncread(momFile,'umean');
V     = ncread(momFile,'vmean');
%
U = mean(U,3,'omitnan');
V = mean(V,3,'omitnan');
ETA=mean(ETA,3,'omitnan');
%
g = 9.8;
[PgrdY,~] = gradientDG(g*ETA./info.dy);
UU = U.*U;
VV = V.*V;
UV = U.*V;
%
H = h+ETA;
%
% get list of velcotiy files
files  = dir([info.rootMat,info.rootName,'uwavg*.nc']);
Nf     = length(files);
%
% preallocate:
uu = 0;
uv = 0;
vv = 0;
uuh= 0;
uvh= 0;
vvh= 0;
for ii=1:Nf
    fprintf('loading from: %s \n', files(ii).name);
    fin = sprintf([info.rootMat,info.rootName,'uwavg_%02d.nc'],ii);
    t   = ncread(fin,'t');
    nt  = length(t);
    %
    % 3) load each variable:
    vars = {'etawavg','uwavg','vwavg'};
    for jj=1:length(vars)
        fin = sprintf([info.rootMat,info.rootName,'%s_%02d.nc'],vars{jj},ii);    
        eval([vars{jj},' = ncread(fin,''',vars{jj},''');'])
    end
    u   = uwavg-U;
    v   = vwavg-V;
    dep = h+etawavg-ETA;
    dep = max(dep,0);
    %
    %
    uu = uu + mean(u.*u,3,'omitnan');
    uv = uv + mean(u.*v,3,'omitnan');
    vv = vv + mean(v.*v,3,'omitnan');
    uuh= uuh+ mean(u.*u.*dep,3,'omitnan');
    uvh= uvh+ mean(u.*v.*dep,3,'omitnan');
    vvh= vvh+ mean(v.*v.*dep,3,'omitnan');
end
uu=uu/Nf;
vv=vv/Nf;
uv=uv/Nf;
uuh=uuh/Nf;
uvh=uvh/Nf;
vvh=vvh/Nf;
%
% spatially smooth all fields (10m in x, 1/2 width of ripchannel in y)
Nflty = 2*info.lc/info.dy; if ~mod(Nflty,2), Nflty=Nflty+1;, end
Nfltx = 10/info.dx;if ~mod(Nfltx,2), Nfltx=Nfltx+1;, end
flt  = hanning(Nflty)*hanning(Nfltx)'; flt = flt./sum(flt(:));
%
uu = conv2(uu,flt,'same');
uv = conv2(uv,flt,'same');
vv = conv2(vv,flt,'same');
%
UU = conv2(UU,flt,'same');
UV = conv2(UV,flt,'same');
VV = conv2(VV,flt,'same');
PgrdY = conv2(PgrdY,flt,'same');
%
% estimate the error in assUming constant depth:
err_uuh = uu-uuh./max(H,1e-1);
err_uvh = uv-uvh./max(H,1e-1);
err_vvh = vv-uvh./max(H,1e-1);
%
% calculate eddy momentum terms
tmp = 0*uu;
tmp(:,2:end-1) = (uu(:,[3:end]) - uu(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uu(:,[2,end]) - uu(:,[1,end-1]))./(1*info.dx);

duudx = tmp;

tmp = 0*uv;
tmp(2:end-1,:) = (uv([3:end],:) - uv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (uv([2,end],:) - uv([1,end-1],:))./(1*info.dy);

duvdy = tmp;

tmp = 0*vv;
tmp(2:end-1,:) = (vv([3:end],:) - vv([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (vv([2,end],:) - vv([1,end-1],:))./(1*info.dy);

dvvdy = tmp;

tmp = 0*uv;
tmp(:,2:end-1) = (uv(:,[3:end]) - uv(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (uv(:,[2,end]) - uv(:,[1,end-1]))./(1*info.dx);

duvdx = tmp;
%
% $$$ tmp = 0*uu;
% $$$ tmp(:,2:end-1) = (err_uuh(:,[3:end]) - err_uuh(:,[1:end-2]))./(2*info.dx);
% $$$ tmp(:,[1 end]) = (err_uuh(:,[2,end]) - err_uuh(:,[1,end-1]))./(1*info.dx);
% $$$ 
% $$$ duudx_err = tmp;
% $$$ 
% $$$ tmp = 0*uv;
% $$$ tmp(2:end-1,:) = (err_uvh([3:end],:) - err_uvh([1:end-2],:))./(2*info.dy);
% $$$ tmp([1 end],:) = (err_uvh([2,end],:) - err_uvh([1,end-1],:))./(1*info.dy);
% $$$ 
% $$$ duvdy_err = tmp;
% $$$ 
% $$$ tmp = 0*vv;
% $$$ tmp(2:end-1,:) = (err_vvh([3:end],:) - err_vvh([1:end-2],:))./(2*info.dy);
% $$$ tmp([1 end],:) = (err_vvh([2,end],:) - err_vvh([1,end-1],:))./(1*info.dy);
% $$$ 
% $$$ dvvdy_err = tmp;
% $$$ 
% $$$ tmp = 0*uv;
% $$$ tmp(:,2:end-1) = (err_uvh(:,[3:end]) - err_uvh(:,[1:end-2]))./(2*info.dx);
% $$$ tmp(:,[1 end]) = (err_uvh(:,[2,end]) - err_uvh(:,[1,end-1]))./(1*info.dx);
% $$$ 
% $$$ duvdx_err = tmp;
%
%
% Calculate mean momentum terms:
tmp = 0*UU;
tmp(:,2:end-1) = (UU(:,[3:end]) - UU(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UU(:,[2,end]) - UU(:,[1,end-1]))./(1*info.dx);

dUUdx = tmp;

tmp = 0*UV;
tmp(2:end-1,:) = (UV([3:end],:) - UV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (UV([2,end],:) - UV([1,end-1],:))./(1*info.dy);

dUVdy = tmp;

tmp = 0*VV;
tmp(2:end-1,:) = (VV([3:end],:) - VV([1:end-2],:))./(2*info.dy);
tmp([1 end],:) = (VV([2,end],:) - VV([1,end-1],:))./(1*info.dy);

dVVdy = tmp;

tmp = 0*UV;
tmp(:,2:end-1) = (UV(:,[3:end]) - UV(:,[1:end-2]))./(2*info.dx);
tmp(:,[1 end]) = (UV(:,[2,end]) - UV(:,[1,end-1]))./(1*info.dx);

dUVdx = tmp;

%
PGY = ncread(momFile,'PgrdY')./H;
RSY = (ncread(momFile,'DySyy') + ncread(momFile,'DxSxy'))./H;
FRY = ncread(momFile,'FRCY')./H;
BRY = ncread(momFile,'BrkDissY')./H;
PGY = conv2(mean(PGY,3,'omitnan'),flt,'same');
RSY = conv2(mean(RSY,3,'omitnan'),flt,'same');
FRY = conv2(mean(FRY,3,'omitnan'),flt,'same');
BRY = conv2(mean(BRY,3,'omitnan'),flt,'same');

