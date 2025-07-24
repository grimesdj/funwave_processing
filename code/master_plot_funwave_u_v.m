%
% get run-directory info
root = '/Users/derekgrimes/funwave/sz2Dturb/';
bath = 'plnr02';
runs  = {'sig04','sig10','sig20','sig40'};
Hs  = 0.5;
Tp  = 8;
%
plotter=1;
%
for ii=1:length(runs)
run = runs{ii};
% get model Hs and surf-zone width
[Hs_x,x,eta_bar,h_bar,sig_eta,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(root,bath,run,Hs,Tp,plotter);
%
% x_range = [xsz+2*(xsz-xsl) xsl];
% now load (u,v,x,y,t)
[u,v,x,y,t] = load_funwave_u_v(root,bath,run);
%
isz = find(x==xsz);
isl = find(x==xsl);
% 
% 
nt = length(t);
nx = length(x);
ny = length(y);
%
% do a little avg
ns = 1*300;
N = floor(nt/ns);
%
f1 = figure;
ha = tight_subplot(1, 2, 0.05, 0.1, 0.1);
axes(ha(1))
ax = gca;
cm = cmocean('balance');
clim1=[-0.05 0.05];
clr1= [clim1(1):diff(clim1)/255:clim1(2)];
%
clim2=[-0.1 0.1];
clr2= [clim2(1):diff(clim2)/255:clim2(2)];
for kk = 1:N
    inds = (kk-1)*ns + (1:ns);
    dum1 = nanmean(u(:,:,inds),3);
    dum2 = nanmean(v(:,:,inds),3);
    axes(ha(1))
    contourf(x,y,dum1,clr1,'edgecolor','none')
    caxis(gca,clim1),colormap(cm)
    ax.Title.String = (sprintf('%02.2f~min',t(inds(1))/60));
    axes(ha(2))
    contourf(x,y,dum2,clr2,'edgecolor','none')
    caxis(gca,clim2),colormap(cm)    
    pause(0.1)
end
%
dum1(dum1>clim(2))=clim(2);
dum1(dum1<clim(1))=clim(1);
dum2(dum2>clim(2))=clim(2);
dum2(dum2<clim(1))=clim(1);
ax.Title.String = sprintf('%02.2f~min',t(inds(1))/60);
axes(ha(1))
contourf(x-xsl,y,dum1,clr,'edgecolor','none')
hold on, plot(xsz*[1 1]-xsl,y([1 end]),'--k','linewidth',2)
%
xlabel('$y~\mathrm{[m]}$','interpreter','latex')
xlabel('$x~\mathrm{[m]}$','interpreter','latex')
set(gca,'color',0.5*[1 1 1],'fontsize',16)
%
a2 = axes;
a2.Position = [ax.Position(1:2)+[0.1 0.4] [0.05 ax.Position(4)/3]];
imagesc(0,clr1,reshape(cm,numel(clr),1,3))
set(a2,'xtick',[])
xlabel('$u$','interpreter','latex','fontsize',18)
%
axes(ha(2))
contourf(x-xsl,y,dum2,clr,'edgecolor','none')
hold on, plot(xsz*[1 1]-xsl,y([1 end]),'--k','linewidth',2)
%
xlabel('$y~\mathrm{[m]}$','interpreter','latex')
% $$$ xlabel('$x~\mathrm{[m]}$','interpreter','latex')
set(gca,'color',0.5*[1 1 1],'fontsize',16,'yticklabel','')
%
a3 = axes;
a3.Position = [ax.Position(1:2)+[0.1 0.4] [0.05 ax.Position(4)/3]];
imagesc(0,clr2,reshape(cm,numel(clr),1,3))
set(a3,'xtick',[])
xlabel('$v$','interpreter','latex','fontsize',18)
%
fout = sprintf('/Users/derekgrimes/funwave/sz2Dturb/mat_data/VORT_%s_%s.mat',bath,run);
save(fout,'omega','x','y','h_bar','xsl','xsz','Lsz')
end
