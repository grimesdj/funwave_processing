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
kp = wavenumber_FunwaveTVD(2*pi*fp,info.hWM);
Ly = 1500;
%
% estimate maximum frequency for hWM
l_min    = pi*info.hWM;
g        = 9.8;
o_max    = sqrt(g*2*pi/l_min*tanh(2*pi/l_min*info.hWM));
f_max    = o_max/(2*pi)
f_valid = find(f<=f_max);
Nfrq_1D = length(f_valid)
%
% first 100 fitting wavenumbers
n_max = min(100,floor(kp*Ly/(2*pi)));
n     = 1:n_max;
%
Ndir = 2*length(n)+1;
Nfrq = floor(nf/Ndir);
df   = range([1/18 f_max])/(Ndir*Nfrq-1);
wave_freq = [1/18:df:f_max]';
%
% $$$ fig1 = figure;
% $$$ fig2 = figure;
% $$$ k = wavenumber_FunwaveTVD(2*pi*wave_freq,info.hWM);
% $$$ theta_n0 = 0*k;
% $$$ for nn = 1:n_max;
% $$$     theta_n = asind( (2*pi*nn/Ly)./k );
% $$$     valid   = imag(theta_n)==0;
% $$$     figure(fig1)
% $$$     hold on,plot(wave_freq(valid),theta_n(valid),'.-',wave_freq(valid),-theta_n(valid),'.-')
% $$$     dtheta   = theta_n-theta_n0;
% $$$     dtheta_theory = (360/Ly)./(k.*cosd(theta_n0));
% $$$     slope = dtheta_theory(valid)\dtheta(valid)
% $$$     theta_n0 = theta_n;
% $$$     figure(fig2)
% $$$     hold on,plot(dtheta_theory(valid),dtheta(valid),'.')
% $$$ end
%
% now begining w/ lowest frequency and using the desired spectral width:
max_theta = 50;
method = 1;
if method == 1
dt0 = 2*max_theta/Ndir;
df0 = (f_max-1/18)/Nfrq;
wave_freq = 1/18;% T = 18s lowest freq
wave_dire = 0;
freq0     = wave_freq(1);
freq_new  = freq0;
theta0    = wave_dire(1);
n = 0;
ind = 1;
while freq_new<f_max
    % estimate the new frequency
    k        = wavenumber_FunwaveTVD(2*pi*freq0,info.hWM);
    dtdn     = (360/Ly)./(k*cosd(max_theta/2));
    for s = [-1 1]
    df1 = df*dtdn/(max_theta/2);
    freq_new = wave_freq(end) + df1;
    wave_dfreq(ind) = df0*dtdn/(max_theta/2);
    % estimate wavenumber
    k_new    = wavenumber_FunwaveTVD(2*pi*freq_new,info.hWM);
    % estimate the approximate new mode number
    if s==-1
        dt = dt0*k_new*Ly*cosd(theta0+dt0)/(360);
        n_new = round(n + dt);
    end
    wave_ddire(ind) = dt;
    % estimate the new direction
    dire_new = asind( s*(2*pi*n_new)/(Ly*k_new) );
    % check if the new theta is outside max_theta, set n=0, and new_freq=wave_freq+delta_f and continue
    if freq_new>f_max
        break
    elseif abs(dire_new)>max_theta | imag(dire_new)~=0
        if s==1
            n = 0;
            theta0 = 0;
            freq0  = freq0 + df0;
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
n_max = min(100,floor(kp*Ly/(2*pi)));
n     = 1:n_max;
%
Ndir = 2*length(n)+1;
Nfrq = floor(nf/Ndir);
df   = range([1/18 f_max])/(Ndir*Nfrq-1);
wave_freq = [1/18:df:f_max]';
%
PHIn = asind(2*pi*n/(Ly*kp)); %atan2( 2*pi*n/Ly, kn)*180/pi;
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
k               = wavenumber_FunwaveTVD(2*pi*wave_freq,info.hWM);
wave_mode       = (Ly/(2*pi)*k.*sind(wave_dire));
wave_mode_up    = ceil(wave_mode);
wave_mode_dwn   = floor(wave_mode);
%
% re-estimate directions
wave_dire_new   = [asin((2*pi*wave_mode_up/Ly)./k)*180/pi, asin((2*pi*wave_mode_dwn/Ly)./k)*180/pi];
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
    return
end
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
    k1 = wavenumber_FunwaveTVD(2*pi*f1,info.hWM);
    for i2 = i1+1:Nf
        f2 = wave_freq(i2);
        d2 = wave_dire(i2);
        o2 = 2*pi*f2;
        a2 = wave_amp(i2);    
        p2 = wave_phi(i2);
        %
        k2 = wavenumber_FunwaveTVD(2*pi*f2,info.hWM);
        %
        f3 = f2-f1;
        d3 = atand( (k2*sind(d2)-k1*sind(d1))/(k2*cosd(d2)-k1*cosd(d1)) );
        cff= cosd(d1-d2+180);
        k3 = sqrt( k1^2 + k2^2 + k1*k2*cff );
        C3 = 2*pi*f3/k3;
        %
        p3 = p1-p2+180;
        p3(p3>=360)=p3(p3>=360)-360;
        %
        % switch sign for estimating coupling coefficient
        o2 = -o2;
        %
        D = (o1+o2).*( o1^2*o2^2/g^2 - k1*k2*cff) - 0.5*( o1*k2^2/cosh(k2*info.hWM) + o2*k1^2/cosh(k1*info.hWM));
        T1 = D*g*(o1 + o2)/(o1*o2)/(g*k3*tanh(k3*info.hWM) - (o1*o2)^2);
        T2 = -g*k1*k2/(2*o1*o2)*cff + 0.5*(o1^2 + o1*o2 + o2^2)/g;
        a3 = abs(T1+T2)*a1*a2;
        %
        nig = nig+1;
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
k_ig = wavenumber_FunwaveTVD(2*pi*ig_freq_band,info.hWM);
k_ig = max(k_ig);
% first 100 fitting wavenumbers
n_max = min(100,floor(k_ig*Ly/(2*pi)));
n     = 1:n_max;
%
Ndir_ig = 2*length(n)+1;
Nfrq_ig = floor(Nf/(2*Ndir));
df_ig   = range(ig_freq_band)/(Ndir_ig*Nfrq_ig-1);
%
dt0_ig = max_theta/Ndir_ig;
df0_ig = diff(ig_freq_band)/Nfrq_ig;
ig_grid_freq = ig_freq_band(1);% T = 18s lowest freq
ig_grid_dire = 0;
freq0_ig     = ig_grid_freq(1);
freq_new  = freq0_ig;
theta0_ig    = ig_grid_dire(1);
n = 0;
ind = 1;
% flag = 0;
while freq_new<ig_freq_band(2);
    % estimate the new frequency
    k        = wavenumber_FunwaveTVD(2*pi*freq_new,info.hWM);
    dtdn     = (360/Ly)./(k);
    for s = [-1 1]
    df1 = df_ig*dtdn/(max_theta/2);
    freq_new = ig_grid_freq(end) + df1;
    ig_grid_dfreq(ind) = df0_ig*dtdn/(max_theta/2);
    % estimate wavenumber
    k_new    = wavenumber_FunwaveTVD(2*pi*freq_new,info.hWM);
    % estimate the approximate new mode number
    if s==-1
        dt = dt0_ig*k_new*Ly*cosd(theta0_ig+dt0_ig)/(360);
        n_new = round(n + max(dt,1));
% $$$         if n == 0 & n_new~=0 & flag == 0
% $$$             n_new = 1;
% $$$             flag  = 1;
% $$$         end
    end
    ig_grid_ddire(ind) = dt;
    % estimate the new direction
    dire_new = asind( s*(2*pi*n_new)/(Ly*k_new) );
    % check if the new theta is outside max_theta, set n=0, and new_freq=ig_freq+delta_f and continue
    if freq_new>ig_freq_band(2)
        break
    elseif abs(dire_new)>max_theta | imag(dire_new)~=0 | dt0_ig/dt>max_theta% | (flag & s==1)
        if s==1
            n = 0;
            %            flag = 0;
            theta0_ig = 0;
            freq0_ig  = freq0_ig + df0_ig;
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
%
%
% next, match ig_vars to nearest ig_grid_vars point and add all amplitudes/phases together.









%
%
% $$$ % first 10 fitting wavenumbers
% $$$ ig_k_avg    = wavenumber_FunwaveTVD(2*pi*ig_freq_avg,info.hWM);
% $$$ n_ig_max       = min(10,floor(ig_k_avg*Ly/(2*pi)));
% $$$ n_ig           = 1:n_max;
% $$$ %
% $$$ Ndir_ig  = 2*n_ig_max+1;
% $$$ Nfrq_ig  = round(nf/(4*Ndir_ig));
% $$$ df_ig    = (min(wave_freq)-min(ig_freq))/(Ndir_ig*Nfrq_ig-1);
% $$$ ig_freq_actual = [min(ig_freq):df_ig:min(wave_freq)]';
% $$$ %
% $$$ PHIn = asind(2*pi*n/(Ly*kp)); %atan2( 2*pi*n/Ly, kn)*180/pi;
% $$$ % each frequency needs to alternate between +/- modes
% $$$ ig_dire_actual = nan(Ndir_ig*Nfrq_ig,1);
% $$$ inds = 1:Ndir_ig*Nfrq_ig;
% $$$ % this is how far each index is from the zero-mode. Even values will be positive, and odd negative.
% $$$ step = round(rem(inds/Ndir_ig,1)*Ndir_ig);
% $$$ zero = step==0;
% $$$ pos  = isodd(step) & ~zero;   
% $$$ neg  = iseven(step)  & ~zero;
% $$$ % give each of these logicals an index corresponding to the mode number
% $$$ ind_pos = (step(pos)+1)/2;
% $$$ ind_neg = (step(neg)/2);
% $$$ ig_dire_actual(zero) = 0;
% $$$ ig_dire_actual(pos)  = PHIn(ind_pos);
% $$$ ig_dire_actual(neg)  =-PHIn(ind_neg);
% $$$ % now estimate wavenumber at each frequency
% $$$ k               = wavenumber_FunwaveTVD(2*pi*ig_freq_actual,info.hWM);
% $$$ ig_mode       = (Ly/(2*pi)*k.*sind(ig_dire_actual));
% $$$ ig_mode_up    = ceil(ig_mode);
% $$$ ig_mode_dwn   = floor(ig_mode);
% $$$ % re-estimate directions
% $$$ ig_dire_new   = [asin((2*pi*ig_mode_up/Ly)./k)*180/pi, asin((2*pi*ig_mode_dwn/Ly)./k)*180/pi];
% $$$ ig_dire_width = diff(ig_dire_new,1,2);
% $$$ % which is closest to original points
% $$$ diff_dire         = abs(ig_dire_actual-ig_dire_new);
% $$$ [diff_dire, best] = min(diff_dire,[],2);
% $$$ ind_best          = sub2ind([Ndir_ig*Nfrq_ig 2],[1:Ndir_ig*Nfrq_ig]',best);
% $$$ ig_dire_actual    = ig_dire_new(ind_best);
%
% first ~10 fitting wavenumbers
n_ig_max = min(10,floor(ig_k_avg*Ly/(2*pi)));
n_ig     = 1:n_ig_max;
%
%
df_ig = range(ig_freq)/(Nfrq_ig-1);
ig_freq_grid = [min(ig_freq):df_ig:max(ig_freq)]';
k_ig         = wavenumber_FunwaveTVD(2*pi*ig_freq_grid,info.hWM);
%
ig_dire_grid  = asin((2*pi*[-n_ig_max:n_ig_max]/Ly)./k_ig)*180/pi;
dd_ig         = gradientDG(ig_dire_grid);
%
%
valid = imag(ig_dire_grid(:))==0;
ig_freq_grid = repmat(ig_freq_grid,Ndir_ig,1);
ig_freq_grid = ig_freq_grid(valid);
ig_dire_grid = ig_dire_grid(valid);
dd_ig        = dd_ig(valid);
N            = sum(valid);
%
%
ig_amp_grid  = nan(N,1);
ig_phi_grid  = nan(N,1);
for ii = 1:N;
        in = find( ig_freq>=ig_freq_grid(ii) - df_ig/2 & ig_freq<ig_freq_grid(ii)     + df_ig/2 &...
                   ig_dire>=ig_dire_grid(ii) - dd_ig(ii)/2 & ig_dire<ig_dire_grid(ii) + dd_ig(ii)/2);
        ig_amp_grid(ii) = sum(ig_amp(in));
        ig_phi_grid(ii) = sum(ig_phi(in));
end
%
%
ig_freq_actual = repmat(ig_freq_grid,Ndir_ig,1).*rand(numel(ig_dire_grid))*df;
valid          = ~isnan(ig_amp_grid) & imag(ig_dire_grid)==0;
ig_freq_actual = ig_freq_grid(valid);
ig_dire_actual = ig_dire_grid(valid);
ig_k_actual    = wavenumber_FunwaveTVD(2*pi*ig_freq_actual',info.hWM);
ig_mode_actual = (Ly/(2*pi)*ig_k_actual.*sind(ig_dire_actual(:)));
ig_dire_actual = asin(round(ig_mode_actual)*(2*pi/Ly)./ig_k_actual )*180/pi;

% add up overlapping (f,d)
[uni_fd,uni] = unique([ig_freq,ig_dire],'rows');
uni_k  =ig_k(uni);
uni_amp=[];
uni_phi=[];
for jj = 1:size(uni_fd,1)
    inds = find( all([ig_freq,ig_dire]==uni_fd(jj,:),2) );
    uni_amp(jj,1) = sum(ig_amp(inds));
    uni_phi(jj,1) = sum(ig_phi(inds));
end
%
% limit amplitudes to >0.5mm
valid_ig = uni_amp>5e-4;
ig_freq = uni_fd(valid_ig,1);
ig_dire = uni_fd(valid_ig,2);
ig_amp  = uni_amp(valid_ig);
ig_phi  = uni_phi(valid_ig);
ig_k    = uni_k(valid_ig);
%
% map to nearest allowable modes
ig_mode       = (Ly/(2*pi)*ig_k.*sind(ig_dire));
ig_mode_up    = ceil(ig_mode);
ig_mode_dwn   = floor(ig_mode);
% re-estimate directions
ig_dire_new   = [asin((2*pi*ig_mode_up/Ly)./ig_k)*180/pi, asin((2*pi*ig_mode_dwn/Ly)./ig_k)*180/pi];
% which is closest to original points
diff_dire     = abs(ig_dire-ig_dire_new);
ig_dire_width = diff(ig_dire_new,1,2); 
[diff_dire, best] = min(diff_dire,[],2);
ind_best          = sub2ind([length(ig_mode) 2],[1:length(ig_mode)]',best);
ig_dire           = ig_dire_new(ind_best);
%
%
valid_ig= abs(ig_dire)<=max_theta & imag(ig_dire)==0;
ig_freq = ig_freq(valid_ig);
ig_dire = real(ig_dire(valid_ig));
ig_amp  = ig_amp(valid_ig);
ig_phi  = ig_phi(valid_ig);
ig_k    = ig_k(valid_ig);
%
