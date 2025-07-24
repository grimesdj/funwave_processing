% decompose depth avg currents into
% rotational and irrotational components
clear all
close all
%
% get run-directory info
root = '/Users/derekgrimes/funwave/sz2Dturb/';
bath = 'plnr02';
runs  = {'sig04','sig10','sig20','sig40'};
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
% x_range = [xsz+2.5*(xsz-xsl) inf];
% now load (u,v,x,y,t)
[u,v,x,y,t] = load_funwave_u_v(root,bath,run,x_range,[-inf inf]);
%
isz = find(x==xsz);
isl = find(x==xsl);
%
dy = median(diff(y));
dx = median(diff(x));
nt = length(t);
%
for jj = 1:nt
    U  = u(:,:,jj);
    V  = v(:,:,jj);
    [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(U,V,dx,dy);
    PSI(:,:,jj) = psi;
    Urot(:,:,jj)= u_psi;
    Vrot(:,:,jj)= v_psi;
    PHI(:,:,jj) = phi;
    Uirr(:,:,jj)= u_phi;
    Virr(:,:,jj)= v_phi;
    %
end
%
fout = sprintf([root,'mat_data/VEL_decomp_%s_%s.mat'],bath,run);
save(fout,'PSI','PHI','Urot','Vrot','Uirr','Virr','x','y','xsl','xsz','Lsz','h_bar','eta_bar','Hs_x','t')
clear PSI PHI Urot Vrot Uirr Virr
end
disp('done!')