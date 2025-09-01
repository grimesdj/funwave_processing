function [out,ky] = alongshore_fft_estimate(info,in);
%
% USAGE: [out,ky] = alongshore_fft_estimate(info,in);
%
% fft of "in" [NxMxP] is over N-rows, M-columns and P-layers are not altered
dy = info.dy;
Ny = info.Ny-1;
%
if ~mod(Ny,2)
    inyq = 1;% don't double the nyquist frequency
    stop = Ny/2+1;
else
    inyq = 0;%
    stop = (Ny+1)/2;
end
ky = [0:stop-1]'/(dy*Ny);
%
%
sz  = size(in);
in  = reshape(in,Ny,prod(sz(2:end)));
out = detrend(in);
out = fft(out);
out = out(1:stop,:);
out = reshape(out,[length(ky), sz(2:end)]);
