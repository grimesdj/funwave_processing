function [Rx,Ry,cFbr,eta,x,y,t,h,mask0] = compile_funwave_breaking_force(rootMat,rootName,bathyFile,subDomain,x,y,h);
%
% USAGE: [Rx,Ry,cFbr,eta,x,y,t,h,mask0] = compile_funwave_breaking_force(rootMat,rootName,bathyFile,subDomain);
%
% load momentum fluxes and breaking eddy viscoscity and estimate breaking wave body force (per unit volume)

% load bathy
%
if ~exist('h','var')
load(bathyFile);
h = -h;% depth needs to be positive!
end
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
%
dx = x(2)-x(1);
dy = y(2)-y(1);
nx = length(x);
ny = length(y);
%
% loop over files and only archive what you need?
files0 = dir([rootMat,rootName,'mask*.mat']);
files1  = dir([rootMat,rootName,'nubrk*.mat']);
files2  = dir([rootMat,rootName,'p*.mat']);
files3  = dir([rootMat,rootName,'q*.mat']);
files4  = dir([rootMat,rootName,'eta*.mat']);
% files4  = dir([rootMat,rootName,'*.mat']);
mask0    = [];
eta0     = [];
nubrk0   = [];
Rx       = [];
Ry       = [];
cFbr     = [];
t0       = [];
for ii=1:length(files0)
    fprintf('loading nubrk from: %s \n', files1(ii).name);
    if exist('subDomain','var')
        dat0 = matfile([files0(ii).folder,'/',files0(ii).name]);
        dat1 = matfile([files1(ii).folder,'/',files1(ii).name]);
        dat2 = matfile([files2(ii).folder,'/',files2(ii).name]);
        dat3 = matfile([files3(ii).folder,'/',files3(ii).name]);
        dat4 = matfile([files4(ii).folder,'/',files4(ii).name]);                        
        mask = dat0.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        nubrk= dat1.nubrk(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        P    = dat2.p(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        Q    = dat3.q(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        eta  = dat4.eta(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);        
        t    = dat1.t;
    else
        load([files0(ii).folder,'/',files0(ii).name]);        
        load([files1(ii).folder,'/',files1(ii).name]);
        load([files2(ii).folder,'/',files2(ii).name]);
        load([files3(ii).folder,'/',files3(ii).name]);
        load([files4(ii).folder,'/',files4(ii).name]);        
    end
    %
    nt = length(t);
    % may need to average to eta-points? Dig into numerics to see where data is output.
    H             = eta + repmat(h,1,1,nt);
    % gradient of fluxes
    [Py, Px, ~]   = gradientDG(P);% (H.*u);
    [Qy, Qx, ~]   = gradientDG(Q);% (H.*v);
    % divergence of fluxes... stress
    [Pyy,   ~, ~] = gradientDG(nubrk.*Py/dy);
    [~  , Pxx, ~] = gradientDG(nubrk.*Px/dx);
    [Qyy,   ~, ~] = gradientDG(nubrk.*Qy/dy);
    [~  , Qxx, ~] = gradientDG(nubrk.*Qx/dx);
    %
    clear Py Px Qy Qx
    Pyy           = Pyy/dy; 
    Pxx           = Pxx/dx;
    Qyy           = Qyy/dy;
    Qxx           = Qxx/dx;
    %
    % body forcing / per unit volume
    Rbx             = (Pyy+Pxx)./H; clear Pxx Pyy
    Rby             = (Qxx+Qyy)./H; clear Qxx Qyy
    [Rbx_y,~    ,~] = gradientDG(Rbx);
    [~    ,Rby_x,~] = gradientDG(Rby);
    %
    % curl of body forcing
    curlFbr = (Rby_x/dx - Rbx_y/dy);
    %
    %
    Rx     = cat(3,Rx,Rbx);
    Ry     = cat(3,Ry,Rby);
    cFbr   = cat(3,cFbr,curlFbr);
    eta0   = cat(3,eta0,eta);
    % new mask estimate using total water depth (h+eta)>0
    %    mask0  = cat(3,mask0,mask);
    mask0  = cat(3,mask0,H>0);    
    nubrk0 = cat(3,nubrk0,nubrk);
    t0     = cat(1,t0,t);
end
eta  = eta0;
mask = mask0;
nubrk = nubrk0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
%
% what should be the "wetted" threshold? (always?, 50%?, other?)
avgmask= mean(mask,3);
mask   = avgmask>0.01;% try 1% wetted!
%
% $$$ % when done processing, use mask to nan bad regions
% $$$ % 5x5 filter
% $$$ NF = 5;
% $$$ filter = ones(NF,NF);filter = filter./sum(filter(:));
% $$$ mask0  = mean(mask,3);
% $$$ mask0  = conv2([ones(2*NF+ny,NF),[ones(NF,nx);mask0;ones(NF,nx)],ones(2*NF+ny,NF)],filter,'same');
% $$$ mask0(mask0<1)=0;
% $$$ mask0 = mask0(NF+1:end-NF,NF+1:end-NF);
%
Rx  = Rx.*mask;
Ry  = Ry.*mask;
cFbr= cFbr.*mask;
eta = eta.*mask;
