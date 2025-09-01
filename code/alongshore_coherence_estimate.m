function [coh,ky,out1,out2,out12] = alongshore_coherence_estimate(info,in1,in2);
dy = info.dy;
Ny = info.Ny;
% need to deal with third dimension... reshape if necessary
sz = size(in1);
if length(sz)<3
    sz(3)=1;
end
in1 = reshape(in1,[sz(1) prod(sz([2 3]))]);
in2 = reshape(in2,[sz(1) prod(sz([2 3]))]);
%
if ~mod(Ny,2)
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
out12= out1.*conj(out2)/(dy*Ny);
out1 = real(out1.*conj(out1))/(dy*Ny);
out2 = real(out2.*conj(out2))/(dy*Ny);
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
%
% account for decreased size of output
sz(1)=stop;
% undo reshape for space/time arrays
coh  = reshape(coh  , sz);
out1 = reshape(out1 , sz);
out2 = reshape(out2 , sz);
out12= reshape(out12, sz);
