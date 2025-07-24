%
% get run-directory info
root = '/Users/derekgrimes/funwave/sz2Dturb/';
bath = 'plnr02';
runs  = {'sig04','sig10','sig20','sig40'};
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
% x_range = [xsz+2.5*(xsz-xsl) inf];
% now load (u,v,x,y,t)
[u,v,x,y,t] = load_funwave_u_v(root,bath,run,x_range,[-inf inf]);
%
isz = find(x==xsz);
isl = find(x==xsl);
% attempt quick curl
[xx,yy] = meshgrid(x,y);
omega = nan(size(u));
% cav   = nan(size(u));
for jj=1:size(u,3);
    [omega(:,:,jj),~] = curl(xx,yy,squeeze(u(:,:,jj)),squeeze(v(:,:,jj)));
end
% 
nt = length(t);
nx = length(x);
ny = length(y);
%
% $$$ % do a 1 min avg
% $$$ ns = 1*60;
% $$$ N = floor(nt/ns);
% $$$ %
% $$$ f1 = figure;
% $$$ ax = axes;
% $$$ cm = cmocean('balance');
% $$$ clim=[-0.1 0.1];
% $$$ clr= [clim(1):diff(clim)/255:clim(2)];
% $$$ for kk = 1:N
% $$$     inds = (kk-1)*ns + (1:ns);
% $$$     dum = nanmean(omega(:,:,inds),3);
% $$$     contourf(x,y,dum,clr,'edgecolor','none')
% $$$     caxis(ax,clim),colormap(ax,cm)
% $$$     ax.Title.String = (sprintf('%02.2f~min',t(inds(1))/60));
% $$$     pause(0.1)
% $$$ end
% $$$ %
% $$$ dum(dum>clim(2))=clim(2);
% $$$ dum(dum<clim(1))=clim(1); 
% $$$ ax.Title.String = sprintf('%02.2f~min',t(inds(1))/60);
% $$$ contourf(x-xsl,y,dum,clr,'edgecolor','none'),
% $$$ caxis(ax,clim),colormap(ax,cm)
% $$$ hold on, plot(xsz*[1 1]-xsl,y([1 end]),'--k','linewidth',2)
% $$$ %
% $$$ xlabel('$y~\mathrm{[m]}$','interpreter','latex')
% $$$ xlabel('$x~\mathrm{[m]}$','interpreter','latex')
% $$$ set(gca,'color',0.5*[1 1 1],'fontsize',16)
% $$$ %
% $$$ a2 = axes;
% $$$ a2.Position = [ax.Position(1:2)+[0.1 0.4] [0.05 ax.Position(4)/3]];
% $$$ imagesc(0,clr,reshape(cm,numel(clr),1,3))
% $$$ set(a2,'xtick',[])
% $$$ xlabel('$\omega_z$','interpreter','latex','fontsize',18)
% $$$ %
% $$$ % save image
% $$$ figname = sprintf([root,bath,'/figures/Vorticity_%s.png'],run);
% $$$ exportgraphics(f1,figname)
% $$$ close(f1)
% save data
fout = sprintf([root,'mat_data/VORT_%s_%s.mat'],bath,run);
save(fout,'omega','x','y','h_bar','xsl','xsz','Lsz')
end
