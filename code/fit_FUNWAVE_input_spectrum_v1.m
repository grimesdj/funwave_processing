function info = fit_FUNWAVE_input_spectrum(info,E,f,d,nf,ig);
%
% USAGE: info = fit_FUNWAVE_input_spectrum(info,nf,ig);
%
% E(f,d): rows are frequency, columns are direction.
% nf: target number of frequency components
% ig: logical for estimating bound IG amplitudes

% we will proceed with the w/ wind version of spectrum, but first need to define the model spectral grid
% get peak frequency
dfrq= f(2)-f(1);% 
ddir= d(2)-d(1);
Ef  = sum(E*ddir,2);
[maxE,imax] = max(Ef);
fp = f(imax);
kp = wavenumber_FunwaveTVD(2*pi*fp,(info.hWM+info.wl));
%
% estimate maximum frequency for hWM
l_min    = pi*(info.hWM+info.wl);
g        = 9.8;
o_max    = sqrt(g*2*pi/l_min*tanh(2*pi/l_min*(info.hWM+info.wl)));
f_max    = o_max/(2*pi);
f_valid = find(f<=f_max);
Nfrq_1D = length(f_valid);
%
% first 100 fitting wavenumbers
n_max = min(100,floor(kp*(info.Ly-info.dy)/(2*pi)));
n     = 1:n_max;
%
Ndir = 2*length(n)+1;
Nfrq = floor(nf/Ndir);
df   = range([1/18 f_max])/(Ndir*Nfrq-1);
wave_freq = [1/18:df:f_max]';
%
% $$$ fig1 = figure;
% $$$ fig2 = figure;
% $$$ k = wavenumber_FunwaveTVD(2*pi*wave_freq,(info.hWM+info.wl));
% $$$ theta_n0 = 0*k;
% $$$ for nn = 1:n_max;
% $$$     theta_n = asind( (2*pi*nn/(info.Ly-info.dy))./k );
% $$$     valid   = imag(theta_n)==0;
% $$$     figure(fig1)
% $$$     hold on,plot(wave_freq(valid),theta_n(valid),'.-',wave_freq(valid),-theta_n(valid),'.-')
% $$$     dtheta   = theta_n-theta_n0;
% $$$     dtheta_theory = (360/(info.Ly-info.dy))./(k.*cosd(theta_n0));
% $$$     slope = dtheta_theory(valid)\dtheta(valid)
% $$$     theta_n0 = theta_n;
% $$$     figure(fig2)
% $$$     hold on,plot(dtheta_theory(valid),dtheta(valid),'.')
% $$$ end
%
% now begining w/ lowest frequency and using the desired spectral width:
max_theta = 50;
f_min = 1/18;
%
method = 1;
if method == 1
dt0 = 2*max_theta/Ndir;
df0 = round(1e5*(f_max-f_min)/Nfrq)*1e-5;
wave_freq = round(1e5*f_min)*1e-5;% T = 18s or ig_band lowest freq
wave_dire = 0;
freq0     = wave_freq(1);
freq_new  = freq0;
theta0    = wave_dire(1);
n = 0;
ind = 1;
while freq_new<f_max
    % estimate the new frequency
    k        = wavenumber_FunwaveTVD(2*pi*freq0,(info.hWM+info.wl));
    dtdn     = (360/(info.Ly-info.dy))./(k*cosd(max_theta/2));
    for s = [-1 1]
    df1 = df*dtdn/(max_theta/2);
    freq_new = round(1e5*(wave_freq(end) + df1))/1e5;
    wave_dfreq(ind) = df0*dtdn/(max_theta/2);
    % estimate wavenumber
    k_new    = wavenumber_FunwaveTVD(2*pi*freq_new,(info.hWM+info.wl));
    % estimate the approximate new mode number
    if s==-1
        dt = dt0*k_new*(info.Ly-info.dy)*cosd(theta0+dt0)/(360);
        n_new = round(n + dt);
    end
    wave_ddire(ind) = dt;
    % estimate the new direction
    dire_new = asind( s*(2*pi*n_new)/((info.Ly-info.dy)*k_new) );
    % check if the new theta is outside max_theta, set n=0, and new_freq=wave_freq+delta_f and continue
    if freq_new>f_max
        break
    elseif abs(dire_new)>max_theta | imag(dire_new)~=0
        if s==1
            n = 0;
            theta0 = 0;
            freq0  = round(1e5*(freq0 + df0))*1e-5;
            ind = ind+1;
            wave_dire(ind,1)= 0;
            wave_freq(ind,1)= freq0;
        end
        continue
    else
        ind = ind+1;
        wave_dire(ind,1) = dire_new;
        wave_freq(ind,1) = freq_new;
        n              = n_new;
    end
    end
end
%
wave_dire_width=wave_ddire';
% $$$ % some points are super close together (where the simplifying assumptions suck)
% $$$ D           = sqrt( (wave_dire-wave_dire').^2/dt0.^2 + (wave_freq-wave_freq').^2/df0.^2);
% $$$ tmp         = 1+diag(nan(size(wave_freq)));
% $$$ [Dmin,imin] = nanmin(D.*tmp,[],2);
% $$$ ibad        = find(Dmin<mean(Dmin)-2*std(Dmin));
% $$$ igood       = [];
% $$$ for ii = 1:length(ibad)
% $$$     isbad = ibad(ii);
% $$$     [pair, ipair] = ismember(imin(isbad),ibad);
% $$$     if pair
% $$$         igood = [igood, ipair];
% $$$     end
% $$$     if ~any(ii==igood)
% $$$         wave_dire(isbad)=nan;
% $$$         wave_freq(isbad)=nan;
% $$$         wave_dfreq(isbad)=nan;
% $$$         wave_ddire(isbad)=nan;
% $$$     end
% $$$ end
%
elseif method==2
%
% first 100 fitting wavenumbers
n_max = min(100,floor(kp*(info.Ly-info.dy)/(2*pi)));
n     = 1:n_max;
%
Ndir = 2*length(n)+1;
Nfrq = floor(nf/Ndir);
df   = range([1/18 f_max])/(Ndir*Nfrq-1);
wave_freq = [1/18:df:f_max]';
%
PHIn = asind(2*pi*n/((info.Ly-info.dy)*kp)); %atan2( 2*pi*n/(info.Ly-info.dy), kn)*180/pi;
% each frequency needs to alternate between +/- modes
wave_dire = nan(Ndir*Nfrq,1);
inds = 1:Ndir*Nfrq;
% this is how far each index is from the zero-mode. Even values will be positive, and odd negative.
step = round(rem(inds/Ndir,1)*Ndir);
zero = step==0;
pos  = isodd(step) & ~zero;   
neg  = iseven(step)  & ~zero;
% give each of these logicals an index corresponding to the mode number
ind_pos = (step(pos)+1)/2;
ind_neg = (step(neg)/2);
wave_dire(zero) = 0;
wave_dire(pos)  = PHIn(ind_pos);
wave_dire(neg)  =-PHIn(ind_neg);
% now estimate wavenumber at each frequency
k               = wavenumber_FunwaveTVD(2*pi*wave_freq,(info.hWM+info.wl));
wave_mode       = ((info.Ly-info.dy)/(2*pi)*k.*sind(wave_dire));
wave_mode_up    = ceil(wave_mode);
wave_mode_dwn   = floor(wave_mode);
%
% re-estimate directions
wave_dire_new   = [asin((2*pi*wave_mode_up/(info.Ly-info.dy))./k)*180/pi, asin((2*pi*wave_mode_dwn/(info.Ly-info.dy))./k)*180/pi];
wave_dire_width = diff(wave_dire_new,1,2);
% which is closest to original points
diff_dire         = abs(wave_dire-wave_dire_new);
[diff_dire, best] = min(diff_dire,[],2);
ind_best        = sub2ind([Ndir*Nfrq 2],[1:Ndir*Nfrq]',best);
wave_dire       = wave_dire_new(ind_best);
wave_dire(zero) = 0;
% populate valid directions
valid          = abs(wave_dire)<=max_theta & imag(wave_dire)==0 & ~isnan(wave_dire);
wave_mode      = wave_mode(valid);
wave_dire      = wave_dire(valid);
wave_freq      = wave_freq(valid);
wave_dire_width= wave_dire_width(valid);
%
% some points are super close together (where the simplifying assumptions suck)
D           = sqrt( (wave_dire-wave_dire').^2/dt0.^2 + (wave_freq-wave_freq').^2/df0.^2);
tmp         = 1+diag(nan(size(wave_freq)));
[Dmin,imin] = nanmin(D.*tmp,[],2);
ibad        = find(Dmin<mean(Dmin)-2*std(Dmin) & wave_dire~=0);
igood       = [];
for ii = 1:length(ibad)
    isbad = ibad(ii);
    [pair, ipair] = ismember(imin(isbad),ibad);
    if pair
        igood = [igood, ipair];
    end
    if ~any(ii==igood)
        wave_dire(isbad)=nan;
        wave_freq(isbad)=nan;
        wave_dire_width(isbad)=nan;
    end
end
%
end
% 
[dd,ff] = meshgrid(d,f);
whos fWM dWM EWM
wave_E = griddata(ff(:)/dfrq,dd(:)/ddir,E(:),wave_freq/dfrq,wave_dire/ddir,'cubic');
wave_E(isnan(wave_E) | wave_E<0)=0;
%
wave_amp = sqrt(wave_E*range(wave_dire)/Ndir*range(wave_freq)/Nfrq.*8)/2;
rng(1)% seeed randum number generator
wave_phi = rand(Nfrq*Ndir,1)*360.0;
%
% remove waves smaller than 1mm
valid    = wave_amp>0.001;
wave_amp = wave_amp(valid);
wave_freq= wave_freq(valid);
wave_dire= wave_dire(valid);
wave_phi = wave_phi(valid);
wave_dire_width = wave_dire_width(valid);
%
%
ig_mode = [];
ig_dire = [];
ig_freq = [];
ig_amp  = [];
ig_phi  = [];
if ~ig
    ig_wave_grid_freq = wave_freq;
    ig_wave_grid_dire = wave_dire;
    ig_wave_grid_amp  = wave_amp;
    ig_wave_grid_phi  = wave_phi;        
else
Nf = length(wave_freq);
%
%
g = 9.8;
nig = 0;
for i1 = 1:Nf-1
    f1 = wave_freq(i1);
    d1 = wave_dire(i1);
    o1 = 2*pi*f1;
    a1 = wave_amp(i1);    
    p1 = wave_phi(i1);
    %
    k1 = wavenumber_FunwaveTVD(2*pi*f1,(info.hWM+info.wl));
    for i2 = i1+1:Nf
        f2 = wave_freq(i2);
        d2 = wave_dire(i2);
        o2 = 2*pi*f2;
        a2 = wave_amp(i2);    
        p2 = wave_phi(i2);
        %
        k2 = wavenumber_FunwaveTVD(2*pi*f2,(info.hWM+info.wl));
        %
        f3 = f2-f1;
        o3 = 2*pi*f3;
        d3 = atand( (k2*sind(d2)-k1*sind(d1))/(k2*cosd(d2)-k1*cosd(d1)) );
        cff= cosd(d1-d2+180);
        k3 = real( sqrt( k1^2 + k2^2 + 2*k1*k2*cff ));
        C3 = 2*pi*f3/k3;
        %
        p3 = p1-p2+180;
        if p3>=360
            p3=p3-360;
        end
        %
        % switch sign for estimating coupling coefficient
        o2 = -o2;
        %
        D = (o1+o2).*( o1^2*o2^2/g^2 - k1*k2*cff) - 0.5*( o1*k2^2/cosh(k2*(info.hWM+info.wl))^2 + o2*k1^2/cosh(k1*(info.hWM+info.wl))^2);
        T1 = D*g*(o1 + o2)/(o1*o2)/(g*k3*tanh(k3*(info.hWM+info.wl)) - (o1+o2)^2);
        T2 = -g*k1*k2/(2*o1*o2)*cff + 0.5*(o1^2 + o1*o2 + o2^2)/g;
        a3 = abs(T1+T2)*a1*a2;
        %
        nig = nig+1;
        if  a3>max(wave_amp)
            disp('issue with coupling coefficient')
            return
        end
        ig_dire(nig,1) = d3;
        ig_freq(nig,1) = f3;
        ig_amp(nig,1)  = a3;
        ig_phi(nig,1)  = p3;
        ig_k(nig,1)    = k3;
    end
end
%
% bin average to nf/4 frequencies and directions
ig_freq_band= [0.005 1/18];
ig_freq_avg = sum(ig_freq_band)/2;
%
% use same spacing as with sea/swell
% k_ig = wavenumber_FunwaveTVD(2*pi*ig_freq_band,(info.hWM+info.wl));
k_ig = kp;% max(k_ig);
% first 100 fitting wavenumbers
n_max = min(100,floor(k_ig*(info.Ly-info.dy)/(2*pi)));
n     = 1:n_max;
%
Ndir_ig = 2*length(n)+1;
Nfrq_ig = floor(Nf/Ndir);
df_ig   = df; % range(ig_freq_band)/(Ndir_ig*Nfrq_ig-1);
%
dt0_ig = dt0; % max_theta/Ndir_ig;
df0_ig = df0; % diff(ig_freq_band)/Nfrq_ig;
ig_grid_freq = round(1e5*(wave_freq(1)-floor((wave_freq(1)-ig_freq_band(1))/df0)*df0))*1e-5;% ig_freq_band(1);% T = 18s lowest freq
ig_grid_dire = 0;
freq0_ig     = ig_grid_freq(1);
freq_new  = freq0_ig;
theta0_ig    = ig_grid_dire(1);
n = 0;
ind = 1;
% flag = 0;
while freq_new<ig_freq_band(2)-df0;
    % estimate the new frequency
    k        = wavenumber_FunwaveTVD(2*pi*freq_new,(info.hWM+info.wl));
    dtdn     = (360/(info.Ly-info.dy))./(k);
    for s = [-1 1]
    df1 = df_ig*dtdn/(max_theta/2);
    freq_new = round(1e5*(ig_grid_freq(end) + df1))*1e-5;
    ig_grid_dfreq(ind) = df0_ig*dtdn/(max_theta/2);
    % estimate wavenumber
    k_new    = wavenumber_FunwaveTVD(2*pi*freq_new,(info.hWM+info.wl));
    % estimate the approximate new mode number
    if s==-1
        dt = dt0_ig*k_new*(info.Ly-info.dy)*cosd(theta0_ig+dt0_ig)/(360);
        n_new = round(n + max(dt,1));
% $$$         if n == 0 & n_new~=0 & flag == 0
% $$$             n_new = 1;
% $$$             flag  = 1;
% $$$         end
    end
    ig_grid_ddire(ind) = dt;
    % estimate the new direction
    dire_new = asind( s*(2*pi*n_new)/((info.Ly-info.dy)*k_new) );
    % check if the new theta is outside max_theta, set n=0, and new_freq=ig_freq+delta_f and continue
    if freq_new>ig_freq_band(2)
        break
    elseif abs(dire_new)>max_theta | imag(dire_new)~=0 | dt0_ig/dt>max_theta% | (flag & s==1)
        if s==1
            n = 0;
            %            flag = 0;
            theta0_ig = 0;
            freq0_ig  = round(1e5*(freq0_ig + df0_ig))*1e-5;
            ind = ind+1;
            ig_grid_dire(ind,1)= 0;
            ig_grid_freq(ind,1)= freq0_ig;
        end
        continue
    else 
        ind = ind+1;
        ig_grid_dire(ind,1) = dire_new;
        ig_grid_freq(ind,1) = freq_new;
        n              = n_new;
    end
    end
end
% combine the ig and wave grids to add ig_amp/phi to grid
N_ig_grid = length(ig_grid_freq);
ig_wave_grid_freq = cat(1,ig_grid_freq,wave_freq);
ig_wave_grid_dire = cat(1,ig_grid_dire,wave_dire);
ig_wave_grid_amp  = cat(1,0*ig_grid_freq,wave_amp);
ig_wave_grid_phi  = cat(1,0*ig_grid_freq,wave_amp);
% next, match ig_vars to nearest ig_grid_vars point and add all amplitudes/phases together.
% should be using dimensionless vars, e.g., ig_freq/df0, etc. 
dT = delaunayTriangulation(ig_wave_grid_freq,ig_wave_grid_dire) ;   % delaunay triangulation 
points = dT.Points ;    % points 
tri = dT.ConnectivityList ;    % nodal connectivity 
x = points(:,1) ;  
X = ig_wave_grid_freq(tri) ;
y = points(:,2) ;
Y = ig_wave_grid_dire(tri) ;
%
%
% then find nearest neighbor grid point, non-dimensionalize to use neighbor distance DI as flag
dT.Points = dT.Points./[df0 dt0];
[NI,DI] = nearestNeighbor(dT,[ig_freq./df0, ig_dire./dt0]);
% unique set of grid points that are neighbors
UI = unique(NI);
i  = sqrt(-1);
% loop over grid points and add all ig_amp/phi
for jj = 1:length(UI)
    ui    = UI(jj);
    iNI   = find(NI==ui & DI<=0.25);
    E0    = ig_wave_grid_amp(ui).*exp(i*ig_wave_grid_phi(ui));
    Eig   = ig_amp(iNI).*exp(i*ig_phi(iNI));
    sumE  = sum(Eig)+E0;
    A     = abs(sumE);
    P     = angle(sumE)*180/pi;
    ig_wave_grid_amp(ui) = A;
    ig_wave_grid_phi(ui) = P;
end
%
% figure, scatter(ig_wave_grid_freq,ig_wave_grid_dire,20,ig_wave_grid_amp,'filled')
% $$$ hold on,plot(wave_freq,wave_dire,'.k','markersize',10)
% $$$        ,plot(ig_grid_freq,ig_grid_dire,'.b','markersize',10)
end
%
%
fig1 = figure;
sc = scatter(ig_wave_grid_dire,ig_wave_grid_freq,18,100*ig_wave_grid_amp(:),'filled','marker','s'); colormap(cmocean('thermal')),set(gca,'color','k','xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
if ig
    figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_scatter_with_boundIG.png'];
else
    figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_scatter.png'];
end
exportgraphics(fig1,figname)
%
TRI = delaunayTriangulation(ig_wave_grid_freq/df0,ig_wave_grid_dire/dt0) ;   % delaunay triangulation 
fig2 = figure;
sc = trisurf(TRI.ConnectivityList,ig_wave_grid_dire,ig_wave_grid_freq,100*ig_wave_grid_amp(:)); colormap(cmocean('thermal')),
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')%,...
                                                                       %   'cameraposition',[0 0.125 35],'cameratarget',[0 0.125 3],'cameraviewangle',0,'cameraupvector',[0 1 0])
view(0,90)
if ig
    figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_with_boundIG.png'];
else
    figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude.png'];
end
exportgraphics(fig2,figname)
%
% $$$ k = wavenumber_FunwaveTVD(2*pi*ig_wave_grid_freq,(info.hWM+info.wl));
% $$$ n = sind(ig_wave_grid_dire).*(info.Ly-info.dy).*k./(2*pi);
% $$$ k1=wavenumber_FunwaveTVD(2*pi*round(1e5*ig_wave_grid_freq)*1e-5,(info.hWM+info.wl));
% $$$ n1= sind(round(1e5*ig_wave_grid_dire)*1e-5).*(info.Ly-info.dy).*k./(2*pi);
% round dire up slightly so that mode is not 0.99999 which maps to 0.0
pm   = sign(ig_wave_grid_dire);
ig_wave_grid_dire = pm.*ceil( 1e5*abs(ig_wave_grid_dire) )*1e-5;
%
%
if ig
    wave_file = [info.rootSim,filesep,'waves_with_boundIG.',info.waveSource];    
else
    wave_file = [info.rootSim,filesep,'waves.',info.waveSource];
end
fid = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(ig_wave_grid_freq));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%2.8f   - Freq \n',ig_wave_grid_freq);
fprintf(fid,'%2.8f   - Dire \n',ig_wave_grid_dire);
fclose(fid);
dlmwrite(wave_file,ig_wave_grid_amp,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,ig_wave_grid_phi,'delimiter','\t','-append','precision',5);
%
% estimate Hs from scattered data
%
info.waveFile = wave_file;
info.waveNdir = length(ig_wave_grid_freq);
info.waveNfrq = length(ig_wave_grid_freq);
info.Tp       = 1/fp;
Hs            = sqrt(8*sum((ig_wave_grid_amp).^2));
info.Hs       = Hs;
info.freqRNG  = [min(wave_freq) max(wave_freq)];
info.direRNG  = [min(wave_dire) max(wave_dire)];
%
save(info.fileName,'-struct','info')
%
%
% now estimate the peak, and spectrum width at each frequncy band then make narrow
% i'm neglecting the IG band here because of the coarse angular resolution
i0 = find(ig_wave_grid_dire==0 & ig_wave_grid_freq>=1/20);
i0 = i0(2:end-1);
i1 = i0+1;
f0 = ig_wave_grid_freq(i0);
n0 = length(f0);
fitfun = fittype( @(a,b,c,x) a*cosd((x-b)/2).^(2*c) );
x0     = [1e-3 0 1];
sigma0 = 10*(pi/180);
sigma  = [];
% tripple the frequency resolution to get roughly 2/3 the same number of wave modes
freq_new    = [ig_wave_grid_freq(1):range(ig_wave_grid_freq)/(3*length(ig_wave_grid_freq)-1):ig_wave_grid_freq(end)]';
narrow_freq = [];
narrow_dire = [];
narrow_amp  = [];
narrow_phi  = [];
% $$$ figure, 
for ii = 1:n0
    f1 = f0(ii);
    i1 = find( ig_wave_grid_freq>=f1 & ig_wave_grid_freq<f1+df0 );
    n1 = length(i1);
    fitout = fit(ig_wave_grid_dire(i1),ig_wave_grid_amp(i1),fitfun,'startpoint',x0);
    %
    sigma(ii,1) = sqrt(2/(fitout.c+1))*(180/pi);
    i2 = find( freq_new>=f1 & freq_new<f1+df0 );
    n2 = length(i2);
% $$$     freq_new = f1+[0:1/(2*n1-1):1]'*df0;
    dire_new = [[2/n2:2/n2:1]*50;[-2/n2:-2/n2:-1]*50];
    dire_new = [0;dire_new(1:n2-1)'];
% $$$     dire_new = fitout.b+acosd( cosd( ig_wave_grid_dire(i1) - fitout.b  ).^(fitout.c/(2/sigma0^2-1) ));
    amp_new  = fitout.a.*cosd( dire_new - fitout.b ).^(2*(2/sigma0^2-1));
% $$$     new = fitout.a.*cosd( ig_wave_grid_dire(i1) - fitout.b ).^(2*(2/sigma0^2-1));
    %
% $$$     plot(ig_wave_grid_dire(i1),ig_wave_grid_amp(i1),'*k')
% $$$     hold on,plot(ig_wave_grid_dire(i1),fitout(ig_wave_grid_dire(i1)),'.r')
% $$$     hold on,plot(ig_wave_grid_dire(i1),new,'.b')
% $$$     pause(0.5)
% $$$     narrow_freq = cat(1,narrow_freq,ig_wave_grid_freq(i1));
% $$$     narrow_dire = cat(1,narrow_dire,ig_wave_grid_dire(i1));
% $$$     narrow_amp  = cat(1,narrow_amp, new);
% $$$     narrow_phi  = cat(1,narrow_phi, ig_wave_grid_phi(i1));    
    narrow_freq = cat(1,narrow_freq,freq_new(i2));    
    narrow_dire = cat(1,narrow_dire,dire_new);
    narrow_amp  = cat(1,narrow_amp, amp_new);
    narrow_phi  = cat(1,narrow_phi, rand(n2,1)*360);
end
%
TRI = delaunayTriangulation(narrow_freq/df0,narrow_dire/dt0) ;   % delaunay triangulation 
fig2 = figure;
sc = trisurf(TRI.ConnectivityList,narrow_dire,narrow_freq,100*narrow_amp(:)); colormap(cmocean('thermal')),
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')%,...
                                                                       %   'cameraposition',[0 0.125 35],'cameratarget',[0 0.125 3],'cameraviewangle',0,'cameraupvector',[0 1 0])
view(0,90)
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_narrow_spectrum.png'];
exportgraphics(fig2,figname)
%
Hs_narrow_spectrum = sqrt(8*sum((narrow_amp).^2))
valid = narrow_amp>0.0001;
narrow_amp = narrow_amp(valid)*Hs./Hs_narrow_spectrum;
narrow_phi = narrow_phi(valid);
narrow_freq = narrow_freq(valid);
narrow_dire = narrow_dire(valid);
%
% try to double frequency resolution for the narrow run
[narrow_freq_double,srt] = sort([narrow_freq; narrow_freq+0.5*gradient(narrow_freq)]);
narrow_dire_double = [narrow_dire; narrow_dire+0.5*gradient(narrow_dire)];
narrow_phi_double  = [narrow_phi ; rand(length(narrow_phi),1)*360];
narrow_dire_double = narrow_dire_double(srt);
narrow_phi_double  = narrow_phi_double(srt);
narrow_amp_double  = griddata(narrow_freq/df0,narrow_dire/dt0,narrow_amp,narrow_freq_double/df0,narrow_dire_double/dt0,'cubic');
%
%
wave_file = [info.rootSim,filesep,'waves_narrow_spectrum.',info.waveSource];
fid = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(narrow_freq));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%8.5f   - Freq \n',narrow_freq);
fprintf(fid,'%10.3f   - Dire \n',narrow_dire);
fclose(fid);
dlmwrite(wave_file,narrow_amp,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,narrow_phi,'delimiter','\t','-append','precision',5);
%
% $$$ %
% $$$ %
% $$$ % $$$ % first 10 fitting wavenumbers
% $$$ % $$$ ig_k_avg    = wavenumber_FunwaveTVD(2*pi*ig_freq_avg,(info.hWM+info.wl));
% $$$ % $$$ n_ig_max       = min(10,floor(ig_k_avg*(info.Ly-info.dy)/(2*pi)));
% $$$ % $$$ n_ig           = 1:n_max;
% $$$ % $$$ %
% $$$ % $$$ Ndir_ig  = 2*n_ig_max+1;
% $$$ % $$$ Nfrq_ig  = round(nf/(4*Ndir_ig));
% $$$ % $$$ df_ig    = (min(wave_freq)-min(ig_freq))/(Ndir_ig*Nfrq_ig-1);
% $$$ % $$$ ig_freq_actual = [min(ig_freq):df_ig:min(wave_freq)]';
% $$$ % $$$ %
% $$$ % $$$ PHIn = asind(2*pi*n/((info.Ly-info.dy)*kp)); %atan2( 2*pi*n/(info.Ly-info.dy), kn)*180/pi;
% $$$ % $$$ % each frequency needs to alternate between +/- modes
% $$$ % $$$ ig_dire_actual = nan(Ndir_ig*Nfrq_ig,1);
% $$$ % $$$ inds = 1:Ndir_ig*Nfrq_ig;
% $$$ % $$$ % this is how far each index is from the zero-mode. Even values will be positive, and odd negative.
% $$$ % $$$ step = round(rem(inds/Ndir_ig,1)*Ndir_ig);
% $$$ % $$$ zero = step==0;
% $$$ % $$$ pos  = isodd(step) & ~zero;   
% $$$ % $$$ neg  = iseven(step)  & ~zero;
% $$$ % $$$ % give each of these logicals an index corresponding to the mode number
% $$$ % $$$ ind_pos = (step(pos)+1)/2;
% $$$ % $$$ ind_neg = (step(neg)/2);
% $$$ % $$$ ig_dire_actual(zero) = 0;
% $$$ % $$$ ig_dire_actual(pos)  = PHIn(ind_pos);
% $$$ % $$$ ig_dire_actual(neg)  =-PHIn(ind_neg);
% $$$ % $$$ % now estimate wavenumber at each frequency
% $$$ % $$$ k               = wavenumber_FunwaveTVD(2*pi*ig_freq_actual,(info.hWM+info.wl));
% $$$ % $$$ ig_mode       = ((info.Ly-info.dy)/(2*pi)*k.*sind(ig_dire_actual));
% $$$ % $$$ ig_mode_up    = ceil(ig_mode);
% $$$ % $$$ ig_mode_dwn   = floor(ig_mode);
% $$$ % $$$ % re-estimate directions
% $$$ % $$$ ig_dire_new   = [asin((2*pi*ig_mode_up/(info.Ly-info.dy))./k)*180/pi, asin((2*pi*ig_mode_dwn/(info.Ly-info.dy))./k)*180/pi];
% $$$ % $$$ ig_dire_width = diff(ig_dire_new,1,2);
% $$$ % $$$ % which is closest to original points
% $$$ % $$$ diff_dire         = abs(ig_dire_actual-ig_dire_new);
% $$$ % $$$ [diff_dire, best] = min(diff_dire,[],2);
% $$$ % $$$ ind_best          = sub2ind([Ndir_ig*Nfrq_ig 2],[1:Ndir_ig*Nfrq_ig]',best);
% $$$ % $$$ ig_dire_actual    = ig_dire_new(ind_best);
% $$$ %
% $$$ % first ~10 fitting wavenumbers
% $$$ n_ig_max = min(10,floor(ig_k_avg*(info.Ly-info.dy)/(2*pi)));
% $$$ n_ig     = 1:n_ig_max;
% $$$ %
% $$$ %
% $$$ df_ig = range(ig_freq)/(Nfrq_ig-1);
% $$$ ig_freq_grid = [min(ig_freq):df_ig:max(ig_freq)]';
% $$$ k_ig         = wavenumber_FunwaveTVD(2*pi*ig_freq_grid,(info.hWM+info.wl));
% $$$ %
% $$$ ig_dire_grid  = asin((2*pi*[-n_ig_max:n_ig_max]/(info.Ly-info.dy))./k_ig)*180/pi;
% $$$ dd_ig         = gradientDG(ig_dire_grid);
% $$$ %
% $$$ %
% $$$ valid = imag(ig_dire_grid(:))==0;
% $$$ ig_freq_grid = repmat(ig_freq_grid,Ndir_ig,1);
% $$$ ig_freq_grid = ig_freq_grid(valid);
% $$$ ig_dire_grid = ig_dire_grid(valid);
% $$$ dd_ig        = dd_ig(valid);
% $$$ N            = sum(valid);
% $$$ %
% $$$ %
% $$$ ig_amp_grid  = nan(N,1);
% $$$ ig_phi_grid  = nan(N,1);
% $$$ for ii = 1:N;
% $$$         in = find( ig_freq>=ig_freq_grid(ii) - df_ig/2 & ig_freq<ig_freq_grid(ii)     + df_ig/2 &...
% $$$                    ig_dire>=ig_dire_grid(ii) - dd_ig(ii)/2 & ig_dire<ig_dire_grid(ii) + dd_ig(ii)/2);
% $$$         ig_amp_grid(ii) = sum(ig_amp(in));
% $$$         ig_phi_grid(ii) = sum(ig_phi(in));
% $$$ end
% $$$ %
% $$$ %
% $$$ ig_freq_actual = repmat(ig_freq_grid,Ndir_ig,1).*rand(numel(ig_dire_grid))*df;
% $$$ valid          = ~isnan(ig_amp_grid) & imag(ig_dire_grid)==0;
% $$$ ig_freq_actual = ig_freq_grid(valid);
% $$$ ig_dire_actual = ig_dire_grid(valid);
% $$$ ig_k_actual    = wavenumber_FunwaveTVD(2*pi*ig_freq_actual',(info.hWM+info.wl));
% $$$ ig_mode_actual = ((info.Ly-info.dy)/(2*pi)*ig_k_actual.*sind(ig_dire_actual(:)));
% $$$ ig_dire_actual = asin(round(ig_mode_actual)*(2*pi/(info.Ly-info.dy))./ig_k_actual )*180/pi;
% $$$ 
% $$$ % add up overlapping (f,d)
% $$$ [uni_fd,uni] = unique([ig_freq,ig_dire],'rows');
% $$$ uni_k  =ig_k(uni);
% $$$ uni_amp=[];
% $$$ uni_phi=[];
% $$$ for jj = 1:size(uni_fd,1)
% $$$     inds = find( all([ig_freq,ig_dire]==uni_fd(jj,:),2) );
% $$$     uni_amp(jj,1) = sum(ig_amp(inds));
% $$$     uni_phi(jj,1) = sum(ig_phi(inds));
% $$$ end
% $$$ %
% $$$ % limit amplitudes to >0.5mm
% $$$ valid_ig = uni_amp>5e-4;
% $$$ ig_freq = uni_fd(valid_ig,1);
% $$$ ig_dire = uni_fd(valid_ig,2);
% $$$ ig_amp  = uni_amp(valid_ig);
% $$$ ig_phi  = uni_phi(valid_ig);
% $$$ ig_k    = uni_k(valid_ig);
% $$$ %
% $$$ % map to nearest allowable modes
% $$$ ig_mode       = ((info.Ly-info.dy)/(2*pi)*ig_k.*sind(ig_dire));
% $$$ ig_mode_up    = ceil(ig_mode);
% $$$ ig_mode_dwn   = floor(ig_mode);
% $$$ % re-estimate directions
% $$$ ig_dire_new   = [asin((2*pi*ig_mode_up/(info.Ly-info.dy))./ig_k)*180/pi, asin((2*pi*ig_mode_dwn/(info.Ly-info.dy))./ig_k)*180/pi];
% $$$ % which is closest to original points
% $$$ diff_dire     = abs(ig_dire-ig_dire_new);
% $$$ ig_dire_width = diff(ig_dire_new,1,2); 
% $$$ [diff_dire, best] = min(diff_dire,[],2);
% $$$ ind_best          = sub2ind([length(ig_mode) 2],[1:length(ig_mode)]',best);
% $$$ ig_dire           = ig_dire_new(ind_best);
% $$$ %
% $$$ %
% $$$ valid_ig= abs(ig_dire)<=max_theta & imag(ig_dire)==0;
% $$$ ig_freq = ig_freq(valid_ig);
% $$$ ig_dire = real(ig_dire(valid_ig));
% $$$ ig_amp  = ig_amp(valid_ig);
% $$$ ig_phi  = ig_phi(valid_ig);
% $$$ ig_k    = ig_k(valid_ig);
% $$$ %
