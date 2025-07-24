function [out,ky] = alongshore_spectra_estimate(info,in);
dy = info.dy;
Ny = info.Ny-1;
if iseven(Ny)
    inyq = 1;% don't double the nyquist frequency
    stop = Ny/2+1;
else
    inyq = 0;%
    stop = (Ny+1)/2;
end
ky = [0:stop-1]'/(dy*Ny);
out = detrend(in);
out = fft(out);
out = out(1:stop,:);
out = real(out.*conj(out));
out(2:end-inyq,:) = 2*out(2:end-inyq,:);

% average over 5-ky bins
Nf  = 5;
ff  = hamming(Nf)./hamming(Nf);
out = conv2(out,ff,'same');
