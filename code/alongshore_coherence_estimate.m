function [coh,ky,out1,out2,out12] = alongshore_coherence_estimate(info,in1,in2);
dy = info.dy;
Ny = info.Ny;
if iseven(Ny)
    inyq = 1;% don't double the nyquist frequency
    stop = Ny/2+1;
else
    inyq = 0;%
    stop = (Ny+1)/2;
end
ky = [0:stop-1]'/(dy*Ny);
out1 = detrend(in1);
out2 = detrend(in2);
out1 = fft(out1);
out2 = fft(out2);
out1 = out1(1:stop,:);
out2 = out2(1:stop,:);
out12= out1.*conj(out2);
out1 = real(out1.*conj(out1));
out2 = real(out2.*conj(out2));
out1(2:end-inyq,:) = 2*out1(2:end-inyq,:);
out2(2:end-inyq,:) = 2*out2(2:end-inyq,:);
out12(2:end-inyq,:) = 2*out12(2:end-inyq,:);
%
coh  = abs(out12)./sqrt(out1.*out2);
%
% average over 5-ky bins
Nf  = 5;
ff  = hamming(Nf); ff = ff./sum(ff);
out1 = conv2(out1,ff,'same');
out2 = conv2(out2,ff,'same');
out12= conv2(out12,ff,'same');