function info = fit_FUNWAVE_input_spectrum(info,E,f,d,nf);
%
% USAGE: info = fit_FUNWAVE_input_spectrum(info,nf);
%
% E(f,d): rows are frequency, columns are direction.
% nf: target number of frequency components

% this script doesn't mess with bound IG modes
% ig=0;
%
% does funwave's wavemaker use the direction waves are going, not coming from?
% d = -d;
%
% estimate E(f) and mean frequency
dfrq= f(2)-f(1);% 
ddir= abs(d(2)-d(1));
Ef  = sum(E*ddir,2);
[maxE,imax] = max(Ef);
fp = f(imax);
% fp = sum(Ef.*f'.*gradient(f))/sum(Ef.*gradient(f));
kp = wavenumber_FunwaveTVD(2*pi*fp,(info.hWM+info.wl));
%
% estimate maximum frequency for k*hWM<=2
k_min    = 2/(info.hWM+info.wl);
g        = 9.8;
o_max    = sqrt(g*k_min*tanh(k_min*(info.hWM+info.wl)));
f_max    = o_max/(2*pi);
%
f_min = 1/18;
I     = find(f>=f_min & f<=f_max);
df    = gradient(f)';
Hs0   = 4*sqrt(sum(Ef(I).*df(I)))
% angular separation between modes at (f_peak,theta=0)
dt   = floor((360/(info.Ly-info.dy))./(kp));
% maximum angle has to be an integer multiple of dt
max_theta = dt*ceil(50/dt);
%
% resulting number of direciton bins (at peak)
Ndir = 2*max_theta/dt;
%
% $$$ %%%%%% TESTING TO INCREASE DENSITY NEAR PEAK AND ALTERNATE MODES AT HIGHER FREQS %%%%%%
% $$$ % $$$ df0 = range([f_min    f_max])/(Ndir*Nfrq-1);
% $$$ % $$$ weight = min(10,(range(Ef)./Ef)); weight = weight./sum(weight);
% $$$ % $$$ df0_stretched = range([f_min, f_max]).*weight;
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % OLD:
% $$$ Nfrq      = floor(nf/Ndir);
% $$$ df        = range([f_min    f_max])/(Ndir*Nfrq-1) ;
% $$$ df0       = df*Ndir;
% $$$ wave_freq =       [f_min:df:f_max]';
% $$$ % round to 1e-5 precission
% $$$ wave_freq = round(1e5*wave_freq)*1e-5;
% $$$ wave_dire = 0*wave_freq;
% $$$ %
% $$$ n_pos    = floor( (Ndir*Nfrq)/2 );
% $$$ n_neg    = floor( (Ndir*Nfrq/2)-1);
% $$$ if iseven(Ndir*Nfrq)
% $$$     n_pos = n_pos-1;
% $$$ end
% $$$ wave_dire(1:2:Ndir*Nfrq)   =     mod([0:dt:(n_pos)*dt],max_theta);
% $$$ wave_dire(2:2:Ndir*Nfrq)   =-(dt+mod([0:dt:(n_neg)*dt],max_theta));
% $$$ %
% $$$ % get nearest mode of each component
% $$$ k         = wavenumber_FunwaveTVD(2*pi*wave_freq,(info.hWM+info.wl));
% $$$ wave_mode = sind(wave_dire)*(info.Ly-info.dy).*k./(2*pi);
% $$$ %
% $$$ % map directions to nearest mode
% $$$ wave_dire = asind( 2*pi*round(wave_mode)./(k*(info.Ly-info.dy)));
% $$$ %
% $$$ % shift slightly above mode because of round-off errors
% $$$ pm = sign(wave_dire);
% $$$ wave_dire = pm.*ceil(1e5*abs(wave_dire))*1e-5;
[wave_freq, wave_dire,df0,Nfrq] = make_FUNWAVE_freq_dire_vectors(nf,(info.Ly-info.dy),(info.hWM+info.wl),f_min,f_max,max_theta,dt);
% 
[dd,ff] = meshgrid(d,f);
% whos fWM dWM EWM
wave_E = griddata(ff(:)/dfrq,dd(:)/ddir,E(:),wave_freq/dfrq,wave_dire/ddir,'cubic');
wave_E(isnan(wave_E) | wave_E<0)=0;
%
wave_amp = sqrt(wave_E*dt*df0.*8)/2;
rng(1)% seed randum number generator
wave_phi = rand(Nfrq*Ndir,1)*360.0;
%
% remove waves smaller than 1mm
valid    = wave_amp>0.0001;
wave_amp_valid = wave_amp(valid);
wave_freq_valid= wave_freq(valid);
wave_dire_valid= wave_dire(valid);
wave_phi_valid = wave_phi(valid);
Hs_full        = sqrt(8*sum((wave_amp_valid).^2))
wave_amp_valid = (Hs0/Hs_full)*wave_amp_valid;
%
%
fig1 = figure;
sc = scatter(wave_dire,wave_freq,18,100*wave_amp(:),'filled','marker','s'); colormap(cmocean('thermal')),set(gca,'color','k','xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_scatter.png'];
exportgraphics(fig1,figname)
%
TRI = delaunayTriangulation(wave_freq/df0,wave_dire/dt) ;   % delaunay triangulation 
fig2 = figure;
sc = trisurf(TRI.ConnectivityList,wave_dire,wave_freq,100*wave_amp(:)); colormap(cmocean('thermal')),
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')%,...
                                                                       %   'cameraposition',[0 0.125 35],'cameratarget',[0 0.125 3],'cameraviewangle',0,'cameraupvector',[0 1 0])
view(0,90)
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude.png'];
exportgraphics(fig2,figname)
%
%
wave_file = [info.rootSim,filesep,'waves.',info.waveSource];
fid       = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(wave_freq_valid));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%2.8f   - Freq \n',wave_freq_valid);
fprintf(fid,'%2.8f   - Dire \n',wave_dire_valid);
fclose(fid);
dlmwrite(wave_file,wave_amp_valid,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,wave_phi_valid,'delimiter','\t','-append','precision',5);
%
% estimate Hs from scattered data
%
info.waveFile = wave_file;
info.waveNdir = length(wave_freq_valid);
info.waveNfrq = length(wave_freq_valid);
info.Tp       = 1/fp;
info.Hs       = Hs0;
info.freqRNG  = [min(wave_freq_valid) max(wave_freq_valid)];
info.direRNG  = [min(wave_dire_valid) max(wave_dire_valid)];
%
save(info.fileName,'-struct','info')
%
%
% Estimate bulk \theta_p and \sigma_\theta from input spectrum
D    = E./Ef;
a1   = sum(D.*cosd(d)*(d(2)-d(1)),2);
b1   = sum(D.*sind(d)*(d(2)-d(1)),2);
theta1 = atand(b1./a1);
a2   = sum(D.*cosd(2*d)*(d(2)-d(1)),2);
b2   = sum(D.*sind(2*d)*(d(2)-d(1)),2);
%
% get averages and estimate mean spread
a1_avg = sum(a1.*Ef)./sum(Ef);
b1_avg = sum(b1.*Ef)./sum(Ef);
theta1_avg = atand(b1_avg./a1_avg);
a2_avg = sum(a2.*Ef)./sum(Ef);
b2_avg = sum(b2.*Ef)./sum(Ef);
%
sigma1     = 180/pi*sqrt( abs(  2*(1 - abs(a1  .* cosd(  theta1    ) + b1 .*  sind(  theta1    )))));
sigma2     = 180/pi*sqrt( abs(0.5*(1 - abs(a2  .* cosd(2*theta1    ) + b2 .*  sind(2*theta1    )))));
sigma1_avg = 180/pi*sqrt( abs(  2*(1 - abs(a1_avg*cosd(  theta1_avg) + b1_avg*sind(  theta1_avg)))));
sigma2_avg = 180/pi*sqrt( abs(0.5*(1 - abs(a2_avg*cosd(2*theta1_avg) + b2_avg*sind(2*theta1_avg)))));
% use a cosine taper to control spread
fitfun = fittype( @(a,b,c,x) a*cosd((x-b)/2).^(2*c) );
x0     = [1e-3 0 1];
%
% get prelim fit
% preallocate output variables
sigma0 = zeros(size(f));
theta0 = zeros(size(f));
E0 = zeros(size(f));
for i = 2:length(f)-1
    dire   = ones(3,1)*d;
    engy   = E([i-1:i+1],:);
    fitout = fit(dire(:),engy(:),fitfun,'startpoint',x0);
    sigma0(i) = sqrt(2/(fitout.c+1))*(180/pi);
    theta0(i) = fitout.b;
    E0(i)     = sum(fitout(dire));
end
%
%
f00     = f(sigma0~=0);
E00     = E0(sigma0~=0);
theta00 = theta0(sigma0~=0);
sigma00 = sigma0(sigma0~=0);
%
%
fig3 = figure;
a1 = subplot(3,1,1);
plot(f00,E00,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,Ef,'.-b','linewidth',1.5,'markersize',8)
ylabel(a1,'$E$ [m$^2$/Hz]','interpreter','latex')
leg = legend('Fit','Obs');
set(leg,'orientation','horizontal')
%
a2 = subplot(3,1,2);
plot(f00,theta00,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,theta1,'.-b','linewidth',1.5,'markersize',8)
yline(sum(theta0.*E0)./sum(E0),'--k')
yline(theta1_avg,'--b')
ylabel(a2,'$\theta_0$ [$^\circ$]','interpreter','latex')
%
a3 = subplot(3,1,3);
plot(f00,sigma00,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,sigma1,'.-b','linewidth',1.5,'markersize',8)
plot(f ,sigma2,'.-r','linewidth',1.5,'markersize',8)
yline(sum(sigma0.*E0)./sum(E0),'--k')
yline(sigma1_avg,'--b')
yline(sigma2_avg,'--r')
ylabel(a3,'$\sigma_\theta$ [$^\circ$]','interpreter','latex')
xlabel(a3,'$f$ [Hz]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max])
set(a2,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max],'ylim',[-25 25],'ytick',-15:15:15)
set(a3,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max],'ylim',[0 50],'ytick',[0:20:40])
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_bulk_statistics_E_theta0_sigma_theta.png'];
exportgraphics(fig3,figname)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Working to efficiently estimate wave_amp, etc.
%%%%   With similar E(f) spectral shape, but specified
%%%%   Spectral width sigma_theta(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) define the target spread values: (201309291200--> observed_sigma_theta ~ 22-25 degrees)
target_sigma = [20 10 5];
%
for iter = 1:length(target_sigma)
% $$$ % 2) for each target_sigma, we must define a set of dire_new,
% $$$ %    while limiting dire <= max_theta *(target_sigma/observed_sigma_theta);
% $$$ new_max_theta = min(max_theta,ceil(4*target_sigma(iter)));
% $$$ % resulting number of direciton bins (at peak)
% $$$ Ndir_new = 2*new_max_theta/dt;
% $$$ % 
% $$$ Nfrq_new  = floor(nf/Ndir_new);
% $$$ df_new    = min(range([f_min, f_max])/(Ndir_new*Nfrq_new-1));
% $$$ df0_new   = df*Ndir_new;
% $$$ freq_new  = [f_min:df_new:f_max]';
% $$$ % round to 1e-5 precission
% $$$ freq_new  = round(1e5*freq_new)*1e-5;
% $$$ dire_new  = 0*freq_new;
% $$$ n_pos    = floor( (Ndir_new*Nfrq_new)/2 );
% $$$ n_neg    = floor( (Ndir_new*Nfrq_new/2)-1);
% $$$ if iseven(Ndir_new*Nfrq_new)
% $$$     n_pos = n_pos-1;
% $$$ end
% $$$ dire_new(1:2:Ndir_new*Nfrq_new)   =     mod([0:dt:(n_pos)*dt],new_max_theta);
% $$$ dire_new(2:2:Ndir_new*Nfrq_new)   =-(dt+mod([0:dt:(n_neg)*dt],new_max_theta));
% $$$ %
% $$$ % get nearest mode of each component
% $$$ k_new     = wavenumber_FunwaveTVD(2*pi*freq_new,(info.hWM+info.wl));
% $$$ wave_mode = sind(dire_new)*(info.Ly-info.dy).*k_new./(2*pi);
% $$$ %
% $$$ % map directions to nearest mode
% $$$ dire_new = asind( 2*pi*round(wave_mode)./(k_new*(info.Ly-info.dy)));
% $$$ %
% $$$ % shift slightly above mode because of round-off errors
% $$$ pm = sign(dire_new);
% $$$ dire_new = pm.*ceil(1e5*abs(dire_new))*1e-5;
%
% now construct D(theta) using cosine taper
sigma = target_sigma(iter)*(pi/180);
Dnew = cosd( d ).^(2*(2/sigma^2-1));
Dnew(abs(d)>75) =0;
Dnew = Dnew./sum(Dnew);
Enew = Ef*Dnew/ddir;
%
max_theta_new = min(max_theta,ceil(2*max_theta*(target_sigma(iter))/sigma2_avg));
[wave_freq, wave_dire,df1,~] = make_FUNWAVE_freq_dire_vectors(nf,(info.Ly-info.dy),(info.hWM+info.wl),f_min,f_max,max_theta_new,dt);
wave_E = griddata(ff(:)/dfrq,dd(:)/ddir,Enew(:),wave_freq/dfrq,wave_dire/ddir,'cubic');
wave_E(isnan(wave_E) | wave_E<=0)=0;
%
wave_amp = sqrt(wave_E*dt*df1.*8)/2;
rng(1)% seed randum number generator
wave_phi = rand(Nfrq*Ndir,1)*360.0;
%
% remove waves smaller than 1mm
valid    = wave_amp>0.0001;
wave_amp_valid = wave_amp(valid);
wave_freq_valid= wave_freq(valid);
wave_dire_valid= wave_dire(valid);
wave_phi_valid = wave_phi(valid);
Hs_spread(iter)= sqrt(8*sum(wave_amp_valid.^2))
wave_amp_valid = (Hs0/Hs_spread(iter))*wave_amp_valid;
%
%
fig1 = figure;
sc = scatter(wave_dire,wave_freq,18,100*wave_amp(:),'filled','marker','s'); colormap(cmocean('thermal')),set(gca,'color','k','xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_sigma',num2str(target_sigma(iter)),'_scatter.png'];
exportgraphics(fig1,figname)
%
TRI = delaunayTriangulation(wave_freq/df0,wave_dire/dt) ;   % delaunay triangulation 
fig2 = figure;
sc = trisurf(TRI.ConnectivityList,wave_dire,wave_freq,100*wave_amp(:)); colormap(cmocean('thermal')),
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')%,...
                                                                       %   'cameraposition',[0 0.125 35],'cameratarget',[0 0.125 3],'cameraviewangle',0,'cameraupvector',[0 1 0])
view(0,90)
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_sigma',num2str(target_sigma(iter)),'.png'];
exportgraphics(fig2,figname)
%
%
wave_file = [info.rootSim,filesep,'waves_sigma',num2str(target_sigma(iter)),'_spectrum.',info.waveSource];
fid       = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(wave_freq_valid));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%2.8f   - Freq \n',wave_freq_valid);
fprintf(fid,'%2.8f   - Dire \n',wave_dire_valid);
fclose(fid);
dlmwrite(wave_file,wave_amp_valid,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,wave_phi_valid,'delimiter','\t','-append','precision',5);
end
return
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %%%%    The Below Code works for estimating spread over bands
% $$$ %%%%    And for interpolating amplitudes to be shorenormal and
% $$$ %%%%    with specified directional spread
% $$$ %%%%    However, it re-defines wave_freq & wave_dire in a weird way.
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use a cosine taper to fit the peak, and spectrum width at each frequncy band to make shorenormal and narrow spectra
% i'm neglecting the IG band here because of the coarse angular resolution
% 1) get the indices where dire=0, 
i0 = find(wave_dire==0 & wave_freq>=f_min);
i0 = i0(2:end-1);
% log these frequencies, and total number of them
f0 = [wave_freq(i0);f_max];
n0 = length(f0);
%
% for each iteration of reducing sigma_theta, increase the number of wave components with theta=0 
f0_new = sort([f0; f0(1:end-1)+diff(f0)/2]);
%
% increase output freq/dire resolution to accomodate narrower spectrum
freq_new = [f0_new + [df/2:df/2:df0/2]]';% [wave_freq(1):range(wave_freq)/(2*length(wave_freq)-1):wave_freq(end)]';
freq_new = freq_new(:);
% build new directions
dire_new = 0*freq_new;
max_theta_new = floor(max_theta/(2*dt))*dt;
n_new = length(freq_new); if ~mod(n_new,2); n_new = n_new-1; end
dire_new(1:2:length(freq_new))   = mod([ 0:dt:floor(n_new*dt/2)],max_theta_new);
dire_new(2:2:length(freq_new))   =-(dt+mod([0:dt:floor(length(freq_new)/2-1)*dt],max_theta_new));
% map directions to nearest mode
k_tmp    = wavenumber_FunwaveTVD(2*pi*freq_new,(info.hWM+info.wl));
n_tmp    = sind(dire_new)*(info.Ly-info.dy).*k_tmp./(2*pi);
dire_new = asind( 2*pi*round(n_tmp)./(k_tmp*(info.Ly-info.dy)));
% shift slightly above mode because of round-off errors
pm = sign(dire_new);
dire_new = pm.*ceil(1e5*abs(dire_new))*1e-5;
%
% preallocate output variables
sigma0 = 0*f0;
theta0 = 0*f0;
E0 = 0*f0;
fitfun = fittype( @(a,b,c,x) a*cosd((x-b)/2).^(2*c) );
x0     = [1e-3 0 1];
%
%
% make a [20 10 5] degree version of these spectra
sig = [20 10 5];
for iter = 1:length(sig)
sigma = sig(iter)*(pi/180);
%
% 
% log the narrow versions
narrow_freq = [];
narrow_dire = [];
narrow_amp  = [];
narrow_phi  = [];
%
% also save a normally incident versions
narrow_normal_amp  = [];
%
if iter==1% if we're just starting, save the full width normally incident
    full_normal_amp = [];
end
% $$$ figure, 
for ii = 1:n0-1
    % get the current frequency band (freq,dire)
    f1 = f0(ii);
    f2 = f0(ii+1);
    i1 = find( wave_freq>=f1 & wave_freq<f2 );
    n1 = length(i1);
    if n1<5, continue, end
    % fit band amplitude
    fitout = fit(wave_dire(i1),wave_amp(i1),fitfun,'startpoint',x0);
    % log the band's spread, mean direction, and energy
    sigma0(ii,1) = sqrt(2/(fitout.c+1))*(180/pi);
    theta0(ii,1) = fitout.b;
    E0(ii,1)     = sum(0.125*(2*wave_amp(i1)).^2);
    % get the higher resolutioon (freq,dire) in current band
    i2 = find( freq_new>=f1 & freq_new<f2 );
    n2 = length(i2);
    % estimate the amplitudes for the narrower spread... can loop here if we want sigma = (5 10 20) degrees.
    amp_new         = fitout.a.*cosd( dire_new(i2) - fitout.b ).^(2*(2/sigma^2-1));
    amp_new_normal  = fitout.a.*cosd( dire_new(i2)            ).^(2*(2/sigma^2-1));
    if iter==1
        amp_new_full_normal =  fitout.a.*cosd(    wave_dire(i1)   ).^(2*fitout.c);
        full_normal_amp     = cat(1,full_normal_amp, amp_new_full_normal);        
    end
    %
% $$$     plot(wave_dire(i1),wave_amp(i1),'*k')
% $$$     hold on,plot(wave_dire(i1),fitout(wave_dire(i1)),'.r')
% $$$     hold on,plot(wave_dire(i1),new,'.b')
% $$$     pause(0.5)
    %    narrow_freq = cat(1,narrow_freq,wave_freq(i1));
    narrow_freq = cat(1,narrow_freq,freq_new(i2));
    narrow_dire = cat(1,narrow_dire,dire_new(i2));
    narrow_amp  = cat(1,narrow_amp, amp_new);
    narrow_normal_amp  = cat(1,narrow_normal_amp, amp_new_normal);
end
narrow_phi = rand(length(narrow_freq),1)*360.0;
%
%
%
Hs_narrow_spectrum = sqrt(8*sum((narrow_amp).^2))
Hs_narrow_normal_spectrum = sqrt(8*sum((narrow_normal_amp).^2))
%
valid        = narrow_amp>0.001;
valid_normal = narrow_normal_amp>0.001;
narrow_amp         = narrow_amp(valid)*Hs_full./Hs_narrow_spectrum;
narrow_normal_amp  = narrow_normal_amp(valid_normal)*Hs_full./Hs_narrow_spectrum;
narrow_normal_phi  = narrow_phi(valid_normal);
narrow_normal_freq = narrow_freq(valid_normal);
narrow_normal_dire = narrow_dire(valid_normal);
narrow_phi = narrow_phi(valid);
narrow_freq = narrow_freq(valid);
narrow_dire = narrow_dire(valid);
%
if iter==1
    Hs_full_normal_spectrum = sqrt(8*sum((full_normal_amp).^2))
    valid_full_normal = full_normal_amp>0.001;
    full_normal_amp  = full_normal_amp(valid_full_normal)*Hs_full./Hs_full_normal_spectrum;
    full_normal_phi  = wave_phi(valid_full_normal);
    full_normal_freq = wave_freq(valid_full_normal);
    full_normal_dire = wave_dire(valid_full_normal);
end
%
%% STOPPED HERE %%%%%
%
%
TRI = delaunayTriangulation(narrow_freq*numel(narrow_freq)/range(narrow_freq),narrow_dire*numel(narrow_dire)/range(narrow_dire)) ;   % delaunay triangulation 
fig2 = figure;
% $$$ sc = trisurf(TRI.ConnectivityList,narrow_dire,narrow_freq,100*narrow_amp(:),'edgecolor','none'); colormap(cmocean('thermal')),
sc = scatter(narrow_dire,narrow_freq,18,100*narrow_amp(:),'filled','marker','s'); colormap(cmocean('thermal')),set(gca,'color','k','xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
view(0,90)
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_sigma',num2str(sig(iter)),'_spectrum.png'];
exportgraphics(fig2,figname)
%
close(fig2)
TRI = delaunayTriangulation(narrow_normal_freq*numel(narrow_normal_freq)/range(narrow_normal_freq),narrow_normal_dire*numel(narrow_normal_dire)/range(narrow_normal_dire)) ;   % delaunay triangulation 
fig2 = figure;
% $$$ sc = trisurf(TRI.ConnectivityList,narrow_normal_dire,narrow_normal_freq,100*narrow_normal_amp(:),'edgecolor','none'); colormap(cmocean('thermal')),
sc = scatter(narrow_normal_dire,narrow_normal_freq,18,100*narrow_normal_amp(:),'filled','marker','s'); colormap(cmocean('thermal')),set(gca,'color','k','xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
ax = gca;
cb = colorbar;
xlabel('$\theta$ [$^\circ$]','interpreter','latex')
ylabel('$f$ [Hz]','interpreter','latex')
ylabel(cb,'$A(\theta,f)$ [cm]','interpreter','latex')
set(cb,'fontsize',10,'ticklabelinterpreter','latex','tickdir','out')
set(ax,'xdir','reverse','ticklabelinterpreter','latex','tickdir','out')
view(0,90)
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_amplitude_sigma',num2str(sig(iter)),'_normal_spectrum.png'];
exportgraphics(fig2,figname)
%
%
f00     = f0(sigma0~=0);
E00     = E0(sigma0~=0);
theta00 = theta0(sigma0~=0);
sigma00 = sigma0(sigma0~=0);
%
%
fig3 = figure;
a1 = subplot(3,1,1);
plot(f00,E00/df0,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,Ef,'.-b','linewidth',1.5,'markersize',8)
ylabel(a1,'$E$ [m$^2$/Hz]','interpreter','latex')
leg = legend('Fit','Obs');
set(leg,'orientation','horizontal')
%
a2 = subplot(3,1,2);
plot(f00,theta00,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,theta1,'.-b','linewidth',1.5,'markersize',8)
yline(sum(theta0.*E0)./sum(E0),'--k')
yline(theta1_avg,'--b')
ylabel(a2,'$\theta_0$ [$^\circ$]','interpreter','latex')
%
a3 = subplot(3,1,3);
plot(f00,sigma00,'.-k','linewidth',2,'markersize',10), hold on
plot(f ,sigma1,'.-b','linewidth',1.5,'markersize',8)
plot(f ,sigma2,'.-r','linewidth',1.5,'markersize',8)
yline(sum(sigma0.*E0)./sum(E0),'--k')
yline(sigma1_avg,'--b')
yline(sigma2_avg,'--r')
ylabel(a3,'$\sigma_\theta$ [$^\circ$]','interpreter','latex')
xlabel(a3,'$f$ [Hz]','interpreter','latex')
set(a1,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max])
set(a2,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max],'ylim',[-25 25],'ytick',-15:15:15)
set(a3,'ticklabelinterpreter','latex','xlim',[f_min 1.1*f_max],'ylim',[0 50],'ytick',[0:20:40])
figname= [info.rootSim,filesep,'figures',filesep,info.rootName,'wavemaker_bulk_statistics_E_theta0_sigma_theta.png'];
exportgraphics(fig3,figname)
%
%
%
%
wave_file = [info.rootSim,filesep,'waves_sigma',num2str(sig(iter)),'_spectrum.',info.waveSource];
fid = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(narrow_freq));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%8.5f   - Freq \n',narrow_freq);
fprintf(fid,'%10.3f   - Dire \n',narrow_dire);
fclose(fid);
dlmwrite(wave_file,narrow_amp,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,narrow_phi,'delimiter','\t','-append','precision',5);
%
wave_file = [info.rootSim,filesep,'waves_sigma',num2str(sig(iter)),'_normal_spectrum.',info.waveSource];
fid = fopen(wave_file,'w');
fprintf(fid,'%5i      - NumFreq  \n',length(narrow_normal_freq));
fprintf(fid,'%10.3f   - PeakPeriod \n',1./fp);
fprintf(fid,'%8.5f   - Freq \n',narrow_normal_freq);
fprintf(fid,'%10.3f   - Dire \n',narrow_normal_dire);
fclose(fid);
dlmwrite(wave_file,narrow_normal_amp,'delimiter','\t','-append','precision',5);
dlmwrite(wave_file,narrow_normal_phi,'delimiter','\t','-append','precision',5);
end
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
