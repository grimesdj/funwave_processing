function info = make_FUNWAVE_input_spectrum_from_SWAN(info)
%
% USAGE: info = make_FUNWAVE_input_spectrum_from_SWAN(info)
%
% extract 2D spectrum from SWAN output files and average over the alongshore extent of the FUNWAVE-TVD wavemaker location. Save WM-spectrum to FUNWAVE I/O files.
%
yyyy = '2013';
time = datenum([yyyy,info.dateTime],'yyyymmddHHMM');% time of simulation
outDIR       = [info.rootSwan, 'output',filesep];
SpecFile     = dir([outDIR,'SP2D.out']);
SpecIterFile = dir([outDIR,'SP2D_iter.out']);
SpecWindFile = dir([outDIR,'SP2D_wind.out']);
Data     = swan_io_spectrum([SpecFile.folder,filesep,SpecFile.name]);
iterData = swan_io_spectrum([SpecIterFile.folder,filesep,SpecIterFile.name]);
windData = swan_io_spectrum([SpecWindFile.folder,filesep,SpecWindFile.name]);
%
E    = Data.VaDens;
Flag = (Data.VaDens==Data.quantity_exception_values);
E(Flag)=nan;
%
iterE    = iterData.VaDens;
iterFlag = (iterData.VaDens==iterData.quantity_exception_values);
iterE(Flag)=nan;
%
windE    = windData.VaDens;
windFlag = (windData.VaDens==windData.quantity_exception_values);
windE(Flag)=nan;
%
x  = Data.x;
y  = Data.y;
f  = Data.frequency;
df = gradient(f)';
d  = 90-Data.directions;
%
Hs     = 4*sqrt(nanmean(nansum(nansum(E,3).*reshape(df,1,length(f),1),2)*5))
iterHs = 4*sqrt(nanmean(nansum(nansum(iterE,3).*reshape(df,1,length(f),1),2)*5))
windHs = 4*sqrt(nanmean(nansum(nansum(windE,3).*reshape(df,1,length(f),1),2)*5))
%
%
Eavg = squeeze(nanmean(E,1));
iterEavg = squeeze(nanmean(iterE,1));
windEavg = squeeze(nanmean(windE,1));
%
%
Ef   = sum(Eavg.*5,2);
Efbc  = sum(iterEavg.*5,2);
% find peaks of Ef and Efa to extract directional distribution at peak
[~,peakE ] = max(Ef);
Ed = Eavg(peakE,:);
[~,peakEbc ] = max(Efbc);
Edbc = iterEavg(peakEbc,:);
% 1) plot E(f,theta) in rectilinear coordinates
fig1 = figure;
ax1  = subplot(3,3,[1 4]);
plot(Efbc,f,'-r',Ef,f,'-k','linewidth',2)
fmin = 1/20;
fmax = 1/4.5;
yline([fmin fmax],'--r')
yline([f(peakE)],':m')
ylabel(ax1,'$f$ [Hz]','interpreter','latex')
xlabel(ax1,'$E(f)$ [m$^2$/Hz]','interpreter','latex')
set(ax1,'xdir','reverse','ticklabelinterpreter','latex','ylim',[0 2*fmax])
%
ax2  = subplot(3,3,[2 3 5 6]);
clrs = [0.001:0.999/255:0.1];
contourf(d,f,Eavg,clrs,'edgecolor','none'),colormap(cmocean('amp')), caxis([0 0.05]), hold on
[C,H]= contour(d',f,Eavg,[0.001 0.01],'k');    
% $$$     imagesc(77.8-dirs',freqs,log10(E')),colormap(cmocean('speed')),caxis([-5 -1])
hold on, 
yline([fmin fmax],'--m')
yline([f(peakE)],':m')
title('SWAN WM-Spectra w/o BC~~~~~~~~~~~~~')
cb = colorbar;
ylabel(cb,'$\mathrm{m^2 / (Hz\,deg)}$','interpreter','latex')
text(-60,0.4,sprintf('w/o BC $$H_\\mathrm{s} = $$ %1.1f m',Hs),'interpreter','latex','fontsize',9)
text(-60,0.35,sprintf('w/  BC $$H_\\mathrm{s} = $$ %1.1f m',iterHs),'interpreter','latex','fontsize',9)
set(ax2,'xlim',[-90 90],'ydir','normal','xticklabel','','yticklabel','','tickdir','out','ytick',get(ax1,'ytick'),'ylim',get(ax1,'ylim'))
%
ax3  = subplot(3,3,[8 9]);
plot(d,Edbc,'-r',d,Ed,'-k','linewidth',2)
xlabel(ax3,'$\theta$ [$^\circ$]','interpreter','latex')
ylabel(ax3,'$E(\theta)$ [m$^2/(^\circ$\,Hz)]','interpreter','latex','rotation',0,'horizontalalignment','right','verticalalignment','cap')
ax3pos = get(ax3,'position');
ax2pos = get(ax2,'position');
ax3pos(3) = ax2pos(3);
hh     = legend('w BC','w/o BC');
set(hh,'edgecolor','none','fontsize',9)                
set(ax3,'tickdir','out','ticklabelinterpreter','latex','xlim',[-90 90],'position',ax3pos)
figname = sprintf([info.rootSim,filesep,'figures',filesep,'swan_output_wavemaker_spectra_zeroBC_%s.png'],datestr(time,30));
if ~exist([info.rootSim,filesep,'figures'],'dir')
    eval(['!mkdir ', [info.rootSim,filesep,'figures']])
end
exportgraphics(fig1,figname)
%
%
%
%
Ef  = sum(iterEavg.*5,2);
Ew  = sum(windEavg.*5,2);
% find peaks of Ef and Efa to extract directional distribution at peak
[~,peakE ] = max(Ef);
Ed = iterEavg(peakE,:);
[~,peakEw ] = max(Ew);
Edw = windEavg(peakEw,:);
% 1) plot E(f,theta) in rectilinear coordinates
fig1 = figure;
ax1  = subplot(3,3,[1 4]);
plot(Ef,f,'-r',Ew,f,'--k','linewidth',2)
fmin = 1/20;
fmax = 1/4.5;
yline([fmin fmax],'--m')
yline([f(peakEw)],':m')
ylabel(ax1,'$f$ [Hz]','interpreter','latex')
xlabel(ax1,'$E(f)$ [m$^2$/Hz]','interpreter','latex')
set(ax1,'xdir','reverse','ticklabelinterpreter','latex','ylim',[0 2*fmax])
%
ax2  = subplot(3,3,[2 3 5 6]);
clrs = [0.001:0.999/255:0.1];
contourf(d,f,windEavg,clrs,'edgecolor','none'),colormap(cmocean('amp')), caxis([0 0.05]), hold on
[C,H]= contour(d',f,iterEavg,[0.001 0.01],'k');    
% $$$     imagesc(77.8-dirs',freqs,log10(E')),colormap(cmocean('speed')),caxis([-5 -1])
hold on, 
yline([fmin fmax],'--m')
yline([f(peakEw)],':m')
title('SWAN WM-Spectra w/ BC~~~~~~~~~~~~~')
cb = colorbar;
ylabel(cb,'$\mathrm{m^2 / (Hz\,deg)}$','interpreter','latex')
text(-70,0.4,sprintf('w/o wind $$H_\\mathrm{s} = $$ %1.1f m',iterHs),'interpreter','latex','fontsize',9)
text(-70,0.35,sprintf('w/  wind $$H_\\mathrm{s} = $$ %1.1f m',windHs),'interpreter','latex','fontsize',9)
set(ax2,'xlim',[-90 90],'ydir','normal','xticklabel','','yticklabel','','tickdir','out','ytick',get(ax1,'ytick'),'ylim',get(ax1,'ylim'))
%
ax3  = subplot(3,3,[8 9]);
plot(d,Ed,'-r',d,Edw,'--k','linewidth',2)
xlabel(ax3,'$\theta$ [$^\circ$]','interpreter','latex')
ylabel(ax3,'$E(\theta)$ [m$^2/(^\circ$\,Hz)]','interpreter','latex','rotation',0,'horizontalalignment','right','verticalalignment','cap')
ax3pos = get(ax3,'position');
ax2pos = get(ax2,'position');
ax3pos(3) = ax2pos(3);
hh     = legend('w/o wind','w/ wind');
set(hh,'edgecolor','none','fontsize',9)                
set(ax3,'tickdir','out','ticklabelinterpreter','latex','xlim',[-90 90],'position',ax3pos)
figname = sprintf([info.rootSim,filesep,'figures',filesep,'swan_output_wavemaker_spectra_withBC_%s.png'],datestr(time,30));
exportgraphics(fig1,figname)
%
%
%
% Continue using the windEavg field
nf   = 1000;% number of sea/swell frequencies
ig   = 1;
E    = windEavg;
% $$$ Ef   = sum(E*(d(2)-d(1)),2);
%
%
info = fit_FUNWAVE_input_spectrum(info,E,f,d,nf);
% info = fit_FUNWAVE_input_spectrum_v1(info,E,f,d,nf,ig);
