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
fin = [root,'mat_data/WaveHeight_',bath,'_',run,'.mat'];
load(fin,'Hs_x','x','eta_bar','h_bar','xsl','xsz','Lsz')
% $$$ [Hs_x,x,eta_bar,h_bar,sig_eta,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(root,bath,run,Hs,Tp,plotter);
%
x_range = [xsz+2*(xsz-xsl) xsl];
iX = find(x>=x_range(1) & x<=x_range(2));
Hs_x = Hs_x(iX);
h_bar = h_bar(iX);
eta_bar = eta_bar(iX);
%
% now load (u,v,x,y,t)
[u,v,x,y,t] = load_funwave_u_v(root,bath,run,x_range);
%
isz = find(x==xsz);
isl = find(x==xsl);
%
% second order SF
n=2;
% avg in time only 
x_range = [];
t_range = [0 inf];
[Sn,r] = calculate_alongshore_structure_function_Sn(u,v,x,y,t,n,t_range,x_range);
%
S0 = (1+6*(ii-1))*2e-4;
S1 = (1+6*(ii-1))*1e-6;
%
f1 = figure('color','w');
ll = loglog(r,(Sn(:,isz)),'ok',r,S0*r.^(2/3),'--r',r,S1*r.^2,'--b','markersize',10,'linewidth',2);
set(ll(1),'markerfacecolor',0.5*[1 1 1])
set(gca,'xlim',[1 r(end)*2])
text(r(1)    ,3*S0  ,'$\Delta_y^{2/3}$','interpreter','latex','color','r','fontsize',18)
text(r(1)*1.3,4*S1  ,'$\Delta_y^{2}$'  ,'interpreter','latex','color','b','fontsize',18)
%
xlabel('$\Delta_y~\mathrm{log}_{10}\mathrm{[m]}$','interpreter','latex')
ylabel('$S_2(\Delta_y)~\mathrm{log}_{10}\mathrm{[m^2\,s^{-2}]}$','interpreter','latex')
%
figname = sprintf([root,bath,'/figures/Surfzone_2ndOrd_StructFun_%s.png'],run);
exportgraphics(f1,figname)
close(f1)
%
fout = sprintf('/Users/derekgrimes/funwave/sz2Dturb/mat_data/S2_%s_%s.mat',bath,run);
save(fout,'Hs_x','x','eta_bar','h_bar','xsl','xsz','Lsz','r','Sn','n')
end
%
% $$$ %
% $$$ %
% $$$ clrs   = cmocean('balance',9);
% $$$ clrs   = clrs([2 3 end-2 end-1],:);
% $$$ fig1 = figure;
% $$$ ax1 = axes;
% $$$ set(gca,'defaultaxescolororder',clrs)
% $$$ loglog(r,Sn(:,isz),'-o',...
% $$$        r,Sn(:,isz+round((isl-isz)/4)),'-d',...
% $$$        r,Sn(:,isz+round((isl-isz)/2)),'-s',...
% $$$        r,Sn(:,isz+round(3*(isl-isz)/4)),'-x')
% $$$ xlabel('$\Delta y~\mathrm{(m)}$','interpreter','latex')
% $$$ ylabel('$S_2~\mathrm{(m/s)^2}$','interpreter','latex')
% $$$ hh = legend('$x=-L_\mathrm{sz}$','$x=-3L_\mathrm{sz}/4$','$x=-L_\mathrm{sz}/2$','$x=-L_\mathrm{sz}/4$')
% $$$ set(hh,'interpreter','latex','location','southeast','fontsize',15)
% $$$ set(gca,'fontsize',16)
