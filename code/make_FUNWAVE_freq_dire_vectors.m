function [wave_freq, wave_dire, dF, Nfrq] = make_FUNWAVE_freq_dire_vectors(nf,width,depth,f_min,f_max,max_theta,dt_min);
%
% USAGE: [wave_freq, wave_dire, dF, dT] = make_FUNWAVE_freq_dire_vectors(nf,width,depth,f_min,f_max,max_theta,dt_min)
%
% estimate the funwave (freq,dire) components with zero coherence following Salatin et al., 2021

% number of direciton bins (at peak)
Ndir = 2*max_theta/dt_min;
% number of frequency bins
Nfrq      = floor(nf/Ndir);
% fresulting requency resolution
df        = range([f_min    f_max])/(Ndir*Nfrq-1) ;
% frequency separation between zero angles
df0       = df*Ndir;
wave_freq =       [f_min:df:f_max]';
% round to 1e-5 precission
wave_freq = round(1e5*wave_freq)*1e-5;
wave_dire = 0*wave_freq;
%
% estimate wavenumber for each frequency
k         = wavenumber_FunwaveTVD(2*pi*wave_freq,depth);
%
n_pos    = floor( (Ndir*Nfrq)/2 );
n_neg    = floor( (Ndir*Nfrq/2)-1);
if iseven(Ndir*Nfrq)
    n_pos = n_pos-1;
end
wave_dire(1:2:Ndir*Nfrq)   =         mod([0:dt_min:(n_pos)*dt_min],max_theta);
wave_dire(2:2:Ndir*Nfrq)   =-(dt_min+mod([0:dt_min:(n_neg)*dt_min],max_theta));
%
% get nearest mode of each component
k         = wavenumber_FunwaveTVD(2*pi*wave_freq,depth);
wave_mode = sind(wave_dire)*width.*k./(2*pi);
%
% map directions to nearest mode
wave_dire = asind( 2*pi*round(wave_mode)./(k*width));
%
% shift slightly above mode because of round-off errors
pm = sign(wave_dire);
wave_dire = pm.*ceil(1e5*abs(wave_dire))*1e-5;
dF = df0;
dT = dt_min;