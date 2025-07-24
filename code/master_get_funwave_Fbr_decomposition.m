% decompose depth avg currents into
% rotational and irrotational components
clear all
close all
%
% get run-directory info
root = '/Users/derekgrimes/funwave/sz2Dturb/';
bath = 'plnr02';
% runs  = {'sig04','sig10','sig20','sig40'};
runs  = {'sig10'};
Hs  = 0.5;
Tp  = 8;
%
plotter=0;
%
for ii=1:length(runs)
run = runs{ii};
% get model Hs and surf-zone width
% $$$ [Hs_x,x,eta_bar,h_bar,sig_eta,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(root,bath,run,Hs,Tp,plotter);
fin = [root,'mat_data/WaveHeight_',bath,'_',run,'.mat'];
load(fin,'Hs_x','x','eta_bar','h_bar','xsl','xsz','Lsz')
%
% limit estimates to surf-zone region
x_range = [-250 10];
iX = find(x-xsl>=x_range(1) & x-xsl<=x_range(2));
eta_bar = eta_bar(iX);
h_bar   = h_bar(iX);
Hs_x    = Hs_x(iX);
%
fin = sprintf([root,'mat_data/Fbreaking_%s_%s.mat'],bath,run);
load(fin,'Rbx','Rby','x','y','h_bar','xsl','xsz','Lsz')
%
[ny,nx,nt] = size(Rbx);
t = 0:1:nt-1;
%
isz = find(x==xsz);
isl = find(x==xsl);
%
dy = median(diff(y));
dx = median(diff(x));
nt = length(t);
%
for jj = 1:nt
    rx  = Rbx(:,:,jj);
    ry  = Rby(:,:,jj);
    [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(rx,ry,dx,dy);
    PSI(:,:,jj) = psi;
    Rxrot(:,:,jj)= u_psi;
    Ryrot(:,:,jj)= v_psi;
    PHI(:,:,jj) = phi;
    Rxirr(:,:,jj)= u_phi;
    Ryirr(:,:,jj)= v_phi;
    %
end
fout = sprintf([root,'mat_data/Fbr_decomp_%s_%s.mat'],bath,run);
save(fout,'PSI','PHI','Rxrot','Ryrot','Rxirr','Ryirr','x','y','xsl','xsz','Lsz','h_bar','eta_bar','Hs_x','t')
clear PSI PHI Rxrot Ryrot Rxirr Ryirr
end
disp('done!')