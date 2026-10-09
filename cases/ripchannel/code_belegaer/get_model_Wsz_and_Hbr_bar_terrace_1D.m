clear all
close all
%
%
runBATHY = 'planar1D';
[HEIGHT,PERIOD,DIR,SPREAD] = ripchannel_parameter_space(runBATHY);
Hs = [];
Tp = [];
for HH=1:length(HEIGHT)
    H = HEIGHT(HH);
    for TT=1:length(PERIOD)
        T = PERIOD(TT);
        for SS=1:length(SPREAD)
            S = SPREAD(SS);               
            for DD=1:length(DIR)
                D = DIR(DD);
                Hs = cat(2,Hs,H);
                Tp = cat(2,Tp,T);                
                %
                % construct identifier
                runWAVES = sprintf('h%02dt%02ds%02dd%02d',H*10,T,S,D);
                %
                % create info structure
                info = ripchannel_run_info(runBATHY,runWAVES);
            end
        end
    end
end
%
info.bathyFile = [info.rootMat,runBATHY,'_depth.mat'];
fin = '/data2/ripchannel/mat_data/planar1D_Hm0.mat';
load(fin);
load(info.bathyFile)
%
% In future, need to normalize the Hm0 by sqrt(2) in the input files
%
% first use linear shoaling to estimate deviation
x = 0.5*(x(1:end-1)+x(2:end));
h = 0.5*(h(1,1:end-1)+h(1,2:end));
xsl = find(h>=0,1,'first');
%
% begin the shoaling at x=400 m
i0 = find(x>400,1,'first');
Hm0= Hm0(:,1:i0);
h  = h(1:i0);
x  = x(1:i0);
%
% $$$ % start with linear group speed
% $$$ om = 2*pi./Tp';
% $$$ k  = wavenumber_FunwaveTVD(om,h);
% $$$ %
% $$$ % need group speed
g = 9.81;
% $$$ a = -0.39;
% $$$ a1 = (a+1/3);
% $$$ %
% $$$ cg = (2*om).^(-1).*( g*h.*(2*k).*(1-a1*(k.*h).^2)./(1-a*(k.*h).^2) + ...
% $$$                     -g*h.*(k).^2.*(2*a1*(k.*h))./(1-a*(k.*h).^2)    + ...
% $$$                      -g*2*a*(h.*k).^3.*(1-a1*(k.*h).^2)./((1-a*(k.*h).^2).^2));
% $$$ %
% $$$ %
% $$$ Hshoal = Hm0(:,end).*sqrt(cg(:,end)./real(cg));
% $$$ Hshoal = Hm0(:,end).*sqrt(sqrt(h(:,end))./sqrt(h));
%
% just use location of peak wave height
[Hbr,imax] = max(Hm0,[],2);
hbr = h(imax)';
gam = Hbr./hbr;
%
%
% find shoreline location
for jj = 1:size(Hm0,1)
    [~,ishore(jj)] = find(Hm0(jj,:)>0.01,1,'first');
end
%
xs = x(ishore);
xbr= x(imax);
%
fprintf('\n the average surfzone width = %1.2f \n',mean(xbr-xs))
%
%
L0  = g.*Tp.^2/(2*pi);
Ir  = info.slope./sqrt(Hbr'./L0);
fprintf('\n the average value of Hbr/hbr = %1.2f \n',mean(gam))
%
figure,
plot(Ir, gam,'.k','markersize',12)
xl = xline(0.4,'--r');
yl = yline(mean(gam),'--b');
leg = legend(yl,sprintf('$\\gamma_\\mathrm{br}=%1.2f$',mean(gam)));
set(leg,'interpreter','latex')
text(0.3, 0.3,'Spilling')
text(0.425, 0.3,'Plunging')
ylabel('$\gamma_\mathrm{br}$ [~]','interpreter','latex')
xlabel('$\xi$ [~]','interpreter','latex')
set(gca,'ylim',[0 0.75],'tickdir','out','ticklabelinterpreter','latex')
figName = [info.rootDAT,filesep,runBATHY,'_breaking_criterion.png'];
exportgraphics(gcf,figName)
%
%
%
% define bar geometr
bar_dep = median(hbr);
fprintf('\naverage depth of breaking is h=%1.2f m\n',mean(hbr));
%
%
frf_bar_amp_ratio = 0.75;
amplitude = bar_dep*frf_bar_amp_ratio;
%
% first, a terraced profile has...
terrace_w   = amplitude/info.slope*exp(-1/2);
%
% next, a bar/trough profile has larger aspect ratio...
bar_amp = 1.25*amplitude;
bar_w   = 0.75*terrace_w;
%
% now find where to put the Gaussian bar
% 1) need initial guess
% $$$ ibar = find(h>bar_dep+amplitude,1,'first');
% $$$ xbar = x(ibar);
% $$$ bar  = bar_amp*exp( -(x-xbar).^2/(2*bar_w^2) );
% $$$ %
% $$$ %
% $$$ iterrace = find(h>bar_dep+amplitude,1,'first');
% $$$ xterrace = x(iterrace);
% $$$ terrace  = amplitude*exp( -(x-xterrace).^2/(2*terrace_w^2) );
% $$$ %
% $$$ % 2) use newton's method to find root using initail guess?
%
% Or iterate of cross-shore location to find where the bar crest is at the desired depth.
planar = -0.03*(x-xsl);
dx     = x(2)-x(1);
nxs    = round(terrace_w/dx);
ixs    = find((x-xsl)>=nxs & (x-xsl)<=x(end)-nxs);
bar_dep_terrace = [];
bar_dep_trough  = [];
for ii = 1:length(ixs)
    % 1) define the current cross-shore location
    xc = x(ixs(ii))-xsl;
    % 2) construct sum of planar w/ gaussian
    zter = planar + amplitude       *exp(-0.5* ((x-(xsl+xc))/terrace_w   ).^2);
    zbar = planar + bar_amp*exp(-0.5* ((x-(xsl+xc))/bar_w).^2);    
    % 3) search for minimum depth within one width of feature
    hter(ixs(ii)) = max(zter( max(ixs(ii)-nxs,1):min(ixs(ii)+nxs,length(x))));
    hbar(ixs(ii)) = max(zbar( max(ixs(ii)-nxs,1):min(ixs(ii)+nxs,length(x))));    
end
idx_ter = find(hter+bar_dep<=0,1,'first');
idx_bar = find(hbar+bar_dep<=0,1,'first');
%
%
loc_terrace = x(idx_ter)-xsl;
loc_bar     = x(idx_bar)-xsl;
loc_fixed   = 125;
fprintf('\nBased on FRF ratio of bar_amp/bar_depth~0.4,\n')
fprintf('setting terrace amplitude to a=%1.2f m\n',amplitude);
fprintf('terrace width: W=%1.2f m\n',terrace_w);
fprintf('terrace crest location is xc= %3.1f m\n',floor(loc_terrace+xsl));

fprintf('\n\n\nbar amplitude: W=%1.2f m\n',bar_amp);
fprintf('bar width: W=%1.2f m\n',bar_w);
fprintf('bar crest location is xc= %3.1f m\n',floor(loc_bar+xsl));
% fprintf('\nGoing to fix this at (xsl-xc)= %3.1f m\n',loc_fixed);
%
%
terrace = amplitude*exp(-0.5* ((x-(xsl+loc_terrace))/terrace_w   ).^2);
bar     = bar_amp  *exp(-0.5* ((x-(xsl+loc_bar))/bar_w).^2);    
%
%
%
figure, plot((x-xsl),-h,'-k',(x-xsl),-(h-bar),'--r',(x-xsl),-(h-terrace),':b','linewidth',2)
xlabel('$(x-x_\mathrm{sl})$ [m]','interpreter','latex')
ylabel('$z$ [m]','interpreter','latex')
yline(-hbr,':k')
set(gca,'xlim',[0 350],'plotboxaspectratio',[1 0.5 1],'tickdir','out',...
        'ytick',[-10:2:0])
figName = [info.rootDAT,filesep,runBATHY,'_bathymetry_with_hbr.png'];
exportgraphics(gcf,figName)
