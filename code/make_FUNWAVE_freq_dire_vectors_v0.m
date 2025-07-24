function [wave_freq, wave_dire, dF, dT] = make_FUNWAVE_freq_dire_vectors(nf,width,depth,f_min,f_max,max_theta,dt_min);
%
% USAGE: [wave_freq, wave_dire, dF, dT] = make_FUNWAVE_freq_dire_vectors(nf,width,depth,f_min,f_max,max_theta,dt_min)
%
df        = range([f_min, f_max])/(nf-1);
wave_freq = [f_min:df:f_max]';
% round to 1e-5 precission
wave_freq = round(1e5*wave_freq)*1e-5;
wave_dire = 0*wave_freq;
%
% estimate wavenumber for each frequency
k         = wavenumber_FunwaveTVD(2*pi*wave_freq,depth);
%
% 1) loop over frequencies, beginning with wave_dire=0
i0        = 1;
f0        = wave_freq(1);
t0        = wave_dire(1);
zero_flag = 0;
mode      = 0;
iter      = 2;
while iter <= nf
% $$$     if iter==32
% $$$         break
% $$$     end
    % see (2.2) below
    if zero_flag
        zero_flag = 0;
        wave_dire(iter)=0;
        dF(i0:iter) = wave_freq(iter)-f0;
        f0   = wave_freq(iter);
        i0   = iter;
        mode = 0;
        iter = iter+1;
        continue
    end
    % 2) increase wave_dire by 1-mode, or #-modes that satisfy dt>=dt_min
    dt   = 0;
    while dt<dt_min
        mode = mode+1;
        dire = asind( 2*pi*mode./(k(iter)*width));
        dt   = dire-abs(wave_dire(iter-1));
    end
    %
    % 2.2) if new mode is > max_theta, go back to wave_dire(i)=0, log dF(i0:i) = wave_freq(i)-wave_freq(i0);    
    if dire>max_theta
        zero_flag=1;
        continue
    end
    % 2.3) log new dire and save dT = wave_dire(i+1)-wave_dire(i)
    wave_dire(iter:iter+1) = [+dire, -dire];
    if wave_dire(iter-1) == 0
        dT(iter-1:iter+1)      = dt;
    else
        dT(iter:iter+1) = dt;
    end
    iter = iter+2;
end
wave_dire = wave_dire(1:nf);