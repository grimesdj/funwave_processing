function info = compare_FUNWAVE_local_vorticity_to_array_estimate(info)
%
% USAGE: info = compare_FUNWAVE_local_vorticity_to_array_estimate(info)
%
% estimate the local vorticity and a 5-m circular array based average vorticity
% and compare to the virtual mooring and observation estimates

% 1) get station/gauge file stats and the depth of the ROD = DEP_ROD_OBS,
load(info.mooringFile,'ROD_MOD','ROD_OBS')
DEP_ROD_OBS = nanmean(ROD_OBS.p);
% 2) locate all cross-shore locations that have H=DEP_ROD_OBS
% 2a) load the grid and averaged waterlevel elevation
load(info.waveStatsFile,'x','y','h','eta_bar')
DEP_MOD = h + eta_bar;
% 2b) find first instance of DEP_MOD>=DEP_ROD_OBS
[xx,yy] = meshgrid(x,y);
[~ ,index]  = min(abs(DEP_MOD-DEP_ROD_OBS),[],2);
IND_ROD_MOD = index;
x_ROD_MOD   = x(index);% this is a bit noisy!
% 2c) get a list of velocity files
filesU = dir([info.rootMat, filesep, info.rootName, 'u*.mat']);
filesV = dir([info.rootMat, filesep, info.rootName, 'v*.mat']);
filesE = dir([info.rootMat, filesep, info.rootName, 'eta*.mat']);
Nf     = length(filesU);
%
% 3) at each cross-shore location, estimate the local voticity, curl( \vec{u} ), and divergence div( \vec{u})
% 3a) loop through each location, load velocity data appropriate for estimating i) local vorticity and ii) area average vorticity
Ny = length(x_ROD_MOD);
%
% ring-of-doom details:
rROD     = 2.5;
L        = 2*rROD*tan(dtheta/2*pi/180);
A        = 0.5*L*rROD*14;
% angular coordinates of each "ADV" on the ROD
theta0   = 0;
dtheta   = 360/14;
thetaADV = theta0 + [0:13]*dtheta;
nx       = ceil(2*rROD/info.dx); if isodd(nx), nx=nx+1; end
ny       = ceil(2*rROD/info.dy); if isodd(ny), ny=ny+1; end
iADV     = (rROD/info.dx)*cosd(thetaADV);
jADV     = (rROD/info.dy)*sind(thetaADV);
% preallocate vars
vort = [];
div  = [];
vort_avg = [];
div_avg  = [];
vort_rod = [];
div_rod  = [];
curlz    = [];
divrg    = [];
for jj = ny/2+1:Ny-ny/2
    % 3b) determine the indices of neighboring points
    y_inds = jj              + [-ny/2:ny/2];
    x_inds = IND_ROD_MOD(jj) + [-nx/2:nx/2];
    DEP    = DEP_MOD(y_inds,x_inds);
    % 3c) loop through files to get all variables
    t = []; u = []; v = []; e = [];
    for kk = 1:Nf
        uin = matfile([filesU(kk).folder,filesep,filesU(kk).name]);
        vin = matfile([filesV(kk).folder,filesep,filesV(kk).name]);
        ein = matfile([filesE(kk).folder,filesep,filesE(kk).name]);
        nt  = length(t);
        t   = cat(1,t, uin.t);
        u   = cat(3,u, uin.u(y_inds,x_inds,:));
        v   = cat(3,v, vin.v(y_inds,x_inds,:));
        e   = cat(3,e, ein.eta(y_inds,x_inds,:));
        Nt  = length(t);
        for tt = nt+1:Nt
            curlz(:,:,tt) = curl(x_inds*info.dx,y_inds*info.dy,u(:,:,tt),v(:,:,tt));
            divrg(:,:,tt) = divergence(x_inds*info.dx,y_inds*info.dy,u(:,:,tt),v(:,:,tt));            
        end
    end
    %
% $$$     % estimate distances to adv locations
% $$$     [iMOD, jMOD] = meshgrid(x_inds-IND_ROD_MOD(jj),y_inds-jj);
% $$$     D = sqrt( (iMOD(:) - iADV).^2 + (jMOD(:) - jADV).^2 );
%
        % estimate distances to adv locations
    [iMOD, jMOD, tMOD] = ndgrid(y_inds-jj,x_inds-IND_ROD_MOD(jj),t-t(1));
    %
    % loop over adv locations and interpolate from grid
    ua = [];
    ur = [];
    for ll = 1:14
        [iOUT,jOUT,tOUT]  = ndgrid(jADV(ll),iADV(ll),t-t(1));
        u1 = interpn(iMOD,jMOD,tMOD,u,iOUT,jOUT,tOUT);
        v1 = interpn(iMOD,jMOD,tMOD,v,iOUT,jOUT,tOUT);
        H1 = interpn(iMOD,jMOD,tMOD,DEP+e,iOUT,jOUT,tOUT);
% $$$         distances   = D(:,ll);
% $$$         if any(distances==0)
% $$$             coloc = find(distances==0);
% $$$             [r,c] = ind2sub(size(xx_inds), coloc);
% $$$             u1    = u(r,c,:);
% $$$             v1    = v(r,c,:);
% $$$             H1    = DEP(r,c)+e(r,c,:);
% $$$         else
% $$$             neighbors   = find(distances<=1);
% $$$             interpWeight= distances(neighbors).^(-2)./sum(distances(neighbors).^(-2));
% $$$             for mm = 1:length(neighbors)
% $$$                 [r,c] = ind2sub(size(xx_inds), neighbors(mm));                
% $$$                 u0(:,mm) = u(r,c,:);
% $$$                 v0(:,mm) = v(r,c,:);
% $$$                 H0(:,mm) = DEP(r,c)+e(r,c,:);
% $$$             end
% $$$             u1    = u0*interpWeight;
% $$$             v1    = v0*interpWeight;
% $$$             H1    = H0*interpWeight;
% $$$             clear u0 v0 H0
% $$$         end
        ua(1:Nt,ll) = u1*cosd(thetaADV(ll)) + v1*sind(thetaADV(ll));
        ur(1:Nt,ll) =-u1*sind(thetaADV(ll)) + v1*cosd(thetaADV(ll));
        H(1:Nt,ll)  = H1;
        clear u1 v1 H1
    end
% $$$     vort_rod(1:Nt,jj) = L/A*sum(ua,2);%;./sum(H,2);
% $$$     div_rod(1:Nt,jj)  = L/A*sum(ur,2);%;./sum(H,2);
    vort_rod(1:Nt,jj) = L/A*sum(ua.*H,2)./sum(H,2);%;./sum(H,2);
    div_rod(1:Nt,jj)  = L/A*sum(ur.*H,2)./sum(H,2);    
    %
    % estimate local and (rectilinear) area averaged
    vort(1:Nt,jj) = curlz(ny/2,nx/2,:);
    div(1:Nt,jj)  = divrg(ny/2,nx/2,:);
    %
    in    = inpolygon(iMOD,jMOD,iADV,jADV);
    H = DEP_MOD(r,c)+e(r,c,:);
    vort_avg(1:Nt,jj) = sum(sum(curlz.*H.*in,1),2)./sum(sum(H.*in,1),2);
    div_avg(1:Nt,jj)  = sum(sum(divrg.*H.*in,1),2)./sum(sum(H.*in,1),2);
    clear H
end
%
% 5) compare them to each other, and compare the y-averaged versions to the observed and station/gauge based estimate
dt = mean(diff(t));
[Sww,f]=welch_method(vort(:,4:end),dt,25,0.5);
[Sww_rod,f]=welch_method(vort_rod(:,4:end),dt,25,0.5);
[Sww_avg,f]=welch_method(vort_avg(:,4:end),dt,25,0.5);
figure, loglog(f,mean(Sww,2),'k',f,mean(Sww_avg,2),'--k',f,mean(Sww_rod,2),'-r',RD(hr).fm,RD(hr).Svort,'-b')