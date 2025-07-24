function [u,v,x,y,h,t,mask0] = load_funwave_u_v(rootMat,rootName,bathyFile,subDomain,x,y,h);
%
% USAGE: [u,v,x,y,h,t,mask0] = load_funwave_u_v(rootMat,rootName,bathyFile,subDomain);
%
%
% load grid and bottom
if ~exist('h','var')
    load(bathyFile);
end
% change convention from h = z_bottom to h=depth
h = -h;
%
if exist('subDomain','var')
    x = x(subDomain(3):subDomain(4));
    y = y(subDomain(1):subDomain(2));
    h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
end
%
% map bottom points to eta points
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
             h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
nx = length(x);
ny = length(y);
dy = median(diff(y));
%
%
% $$$ if y_range(2)==ny
% $$$     % get nearest power of 2 in alongshore 
% $$$     N2 = log2(y(end)-y(1));
% $$$     i2 = (2^floor(N2)/dy);
% $$$     iY = [1:i2]';
% $$$     y  = y(iY);
% $$$ else
% $$$ iY = find(y>=y_range(1) & y<=y_range(2));
% $$$ Y  = y(iY);
% $$$ end
% $$$ %
% $$$ iX = find(x>=x_range(1) & x<=x_range(2));
% $$$ x  = x(iX);
%
% load data
% $$$ files = dir([rootMat,rootName,'mask*.mat']);
files1= dir([rootMat,rootName,'eta*.mat']);
% $$$ mask0  = [];
eta0  = [];
t0    = [];
for ii=1:length(files1)
    fprintf('loading eta for masking from: %s \n', files1(ii).name);
    if exist('subDomain','var')
% $$$         dat = matfile([files(ii).folder,'/',files(ii).name]);
        dat1= matfile([files1(ii).folder,'/',files1(ii).name]);
% $$$         mask = dat.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        eta  = dat1.eta(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        t    = dat1.t;
    else
% $$$         load([files(ii).folder,'/',files(ii).name]);
        load([files1(ii).folder,'/',files1(ii).name]);        
    end
% $$$     mask0 = cat(3,mask0,mask);
    eta0 = cat(3,eta0,eta);        
    t0   = cat(1,t0,t);
end
% mask = mask0;
eta  = eta0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
%
% make a mask for (u,v) based on wetted region
mask0  = (h + eta>0);
% what should be the "wetted" threshold? (always?, 50%?, other?)
avgmask= mean(mask0,3);
mask   = avgmask>0.01;% try 1% wetted!
% $$$ % 5x5 filter for estimating smooth masking region
% $$$ NF = 5;
% $$$ filter = ones(NF,NF);filter = filter./sum(filter(:));
% $$$ mask0  = mean(mask,3);
% $$$ mask0  = conv2([ones(2*NF+ny,NF),[ones(NF,nx);mask0;ones(NF,nx)],ones(2*NF+ny,NF)],filter,'same');
% $$$ mask0(mask0<1)=0;
% $$$ mask0 = mask0(NF+1:end-NF,NF+1:end-NF);
%
% load u data
files = dir([rootMat,rootName,'u_*.mat']);
u0  = [];
for ii=1:length(files)
    fprintf('loading u from: %s \n', files(ii).name);
    if exist('subDomain','var')
        dat = matfile([files(ii).folder,'/',files(ii).name]);
        u = dat.u(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
    else
        load([files(ii).folder,'/',files(ii).name],'u');
    end
    u0 = cat(3,u0,u);
end
u = u0.*mask;
%
% load v data
files = dir([rootMat,rootName,'v_*.mat']);
v0  = [];
for ii=1:length(files)
    fprintf('loading v from: %s \n', files(ii).name);        
    if exist('subDomain','var')
        dat = matfile([files(ii).folder,'/',files(ii).name]);
        v = dat.v(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
    else
        load([files(ii).folder,'/',files(ii).name],'v');
    end
    v0 = cat(3,v0,v);
end
v = v0.*mask;
%
% $$$ u(~mask)=nan;
% $$$ v(~mask)=nan;