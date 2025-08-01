function [mean_stats,binned_stats,front_log] = compile_funwave_bore_fronts(rootMat,rootName,bathyFile,subDomain,x,y,h,Xbins);
%
% USAGE: [mean_stats,binned_stats,front_log] = compile_funwave_bore_fronts(rootMat,rootName,bathyFile,subDomain,x,y,h,Xbins);
%
% OLD: [mL,sL,N,mLx,sLx,Lx_log_mean,Lx_log_std,Nx,xylog,Xbins] = compile_funwave_bore_fronts(rootMat,rootName,bathyFile,subDomain,x,y,h,Xbins);
%
%
% load grid and bottom
if ~exist('h','var')
   load(bathyFile);
   h = -h;
end
%
if exist('subDomain','var')
    x = x(subDomain(3):subDomain(4));
    y = y(subDomain(1):subDomain(2));
    h = h(subDomain(1):subDomain(2),subDomain(3):subDomain(4));
end
%
%
% map bottom points to eta points
h  = 0.25*(h(1:end-1,1:end-1) + h(2:end,1:end-1) + ...
           h(2:end,1:end-1) + h(2:end,2:end));
x  = 0.5*(x(1:end-1) + x(2:end));
y  = 0.5*(y(1:end-1) + y(2:end));
% make (x,y) column vectors
x  = x(:);
y  = y(:);
nx = length(x);
ny = length(y);
%
%
[yy,xx]   = meshgrid(y',x');
%
dx = x(2)-x(1);
%
% load data
files0 = dir([rootMat,rootName,'mask*.mat']);
files = dir([rootMat,rootName,'nubrk*.mat']);
mask0  = [];
nubrk0   = [];
t0     = [];
for ii=1:length(files)
    fprintf('loading nubrk from: %s \n', files(ii).name);
    if exist('subDomain','var')
        dat0= matfile([files0(ii).folder,'/',files0(ii).name]);
        dat = matfile([files(ii).folder,'/',files(ii).name]);
        mask= dat0.mask(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        nubrk = dat.nubrk(subDomain(1):subDomain(2)-1,subDomain(3):subDomain(4)-1,:);
        t = dat.t;
    else
        load([files0(ii).folder,'/',files0(ii).name]);        
        load([files(ii).folder,'/',files(ii).name]);
    end
    mask0 = cat(3,mask0,mask);    
    nubrk0 = cat(3,nubrk0,nubrk);
    t0   = cat(1,t0,t);
end
mask = mask0;
nubrk = nubrk0;
t   = t0;
nt  = length(t);
dt  = mean(diff(t));
% small search radius as the signal is binary
r0 = floor(2.5/dx);% this is 2.5m in x, and 5m in y for (dx=0.5,dy=1) meters
%
% for cross shore bin averages
if ~exist('Xbins','var');
    db = 2;
    Xbins = [25:db:250];
else
    db = Xbins(2)-Xbins(1);
end
iter  = 1;
for kk = 1:length(t);
    nu0  = nubrk(:,:,kk);
    nu   = nu0~=0;
    rclog= bore_front_search_funwave_v2(nu',ny,nx,r0,0.5);
        if isempty(rclog)
            continue
	end
        for ww = 1:length(rclog)
            % convert row/col to x/y
            rc = rclog{ww};
            if isempty(rc)
                continue
            end
            cf = rc(:,2);
	    rf = rc(:,1);
            xf = x(rf);
            yf = y(cf);
            % make sure points are oriented continuously south to north
            [Y,srt]=sort(yf);
            X = xf(srt);
            % estimate front length
% $$$             dl = sqrt(diff(X).^2 + diff(Y).^2);
% $$$             l(iter)  = sum(dl);
            l(iter)  = max(Y)-min(Y);
            xl(iter) = mean(X);
            yl(iter) = mean(Y);
            nl(iter) = length(X);
            tl(iter) = t(kk);
            iter = iter+1;
            xylog{ww,kk} = [X Y];
% $$$     [xylog,bblog,fnlog] = bore_front_search_v5(Q0,nx,ny,r0,Qmin,Qstd,Nstd,plotter,ax);%        rclog= bore_front_search_v2(nu',ny,nx,r0,0,0.1);
% $$$ % $$$         imagesc(x,y,nu)
% $$$     if isempty(crlog)
% $$$         continue
% $$$     end
% $$$     %
% $$$     for tt = 1:Nt
% $$$         if isempty(crlog{tt}),
% $$$             disp(['no fronts for tp = ',num2str(tt),' in file: ',fin])
% $$$             continue
% $$$         end
% $$$         C = crlog{tt}(:,1);
% $$$         R = crlog{tt}(:,2);
% $$$         ind = sub2ind(size(xx),R, C);
% $$$         X = xx(ind);
% $$$         Y = yy(ind);
% $$$         %
% $$$         % make sure points are oriented continuously south to north 
% $$$         [Y,srt]=sort(Y);
% $$$         X = X(srt);
% $$$         % estimate front length
% $$$         dl = sqrt(diff(X).^2 + diff(Y).^2);
% $$$         L(ww) = sum(dl);
% $$$         % recall that xP is transposed relative to x; 
% $$$         rcF{tt} = [C R];
% $$$         xyF{tt} = [X Y];
% $$$         tF      = [tt t(tt)]+0*R;
    end
end
%
% $$$ front_log = struct('Xavg',xl,'Yavg',yl,'Length',l,'NumPoints',nl,'time',tl,'xylog',xylog);
%
% keep stats
nframes = length(t);
N       = size(xylog,1)/nframes;
mL      = exp(mean(log(l)));
sL      = exp(mean(log(l))+std(log(l)));
Lbins   = [0:10:range(y)];
% histogram of lengths
pL      = hist(l,Lbins);
mean_stats = struct('N',N,'Length',mL,'Length_plus_std',sL,'Lbins',Lbins,'Length_histogram',pL);
%
% cross-shore bins
Nx  = nan*Xbins;
mLx = nan*Xbins;
sLx = nan*Xbins;
Lx_log_mean = nan*Xbins;
Lx_log_std  = nan*Xbins;
Lx_histo    = nan*(Lbins'*Xbins);
for bin = 1:length(Xbins)
    iX = find(xl>=Xbins(bin)-db/2 & xl<Xbins(bin)+db/2);
    Nx(bin) = length(iX)/nframes;
    Lx_log_mean(bin)=mean(log(l(iX)));
    Lx_log_std(bin)=std(log(l(iX)));
    mLx(bin)=exp(mean(log(l(iX))));
    sLx(bin)=exp(log(mLx(bin))+std(log(l(iX))));
    tmp = hist(l(iX),Lbins);
    Lx_histo(:,bin) = tmp;
end
%
%
binned_stats = struct('Xbins',Xbins,'N',Nx,'Length',mLx,'Length_plus_std',sLx,'Length_histogram',Lx_histo,'log_mean_length',Lx_log_mean,'log_std_length',Lx_log_std);

