%
if ~exist('runDay','var')
    root   = '/home/derek/projects/ShortCrests/MOD/funwave/';
    runDay = '0929';
    runID  = 'run00';
end
%
fprintf('processing case: %s  %s \n',runDay, runID)
[rootSim,rootOut,rootMat,rootName,bathyFile,t0,dt0,x0,dx0,y0,dy0,h0,Tp0,Hs0,Dir,Sig,Gam,xWM,hWM] = funwave_run_info(runDay,runID);
[rootSim,rootOut,rootMat,rootName,bathyFile,t,dt,x,dx,y,dy,h,Tp0,Hs0,Dir,Sig,Gam] = funwave_run_info(runDay,runID);
%
% only analyze obs along the 50's line of instruments
iY = find(abs(y-730)<=10);
nx = length(x0);
ny = length(iY);
SubDomain = [iY(1) iY(1)+ny 1 nx];
y  = y0(iY);
h  = h0(iY,:);
plotter=1;
%
[Hs_xy,eta_bar,mask0,x,y,h,freqs,Snn_xy,xsl] = calculate_funwave_wave_statistics(rootMat,rootName,bathyFile,Hs0,Tp0,subDomain,x0,y0,h0);
Hs_x = nanmean(Hs_xy,1);
Snn_x= nanmean(Snn_xy,1);
% $$$ [Hs_x,x,eta_bar,h_bar,sig_eta,freqs,Snn_xy,Snn_x,Snn_wg,xsl,xsz,Lsz] = calculate_funwave_wave_height_statistics(rootMat,rootName,bathyFile,Hs0,Tp0);
%
% use linear theory to shoal spectrum
g = 9.8;
om = 2*pi*freqs;
k0 = wavenumber(om',-h_bar);
df = freqs(2)-freqs(1);
dpth=ones(length(om),1)*max((eta_bar-h_bar),0.1);
Cg = 0.5*(g*tanh(k0.*dpth)+g*k0.*dpth.*(sech(k0.*dpth).^2))./sqrt(g*k0.*tanh(k0.*dpth));
[~,iwg] = min(abs(h_bar-0.99*min(h_bar)));
Snn_linear = (Snn_x(:,iwg).*Cg(:,iwg).*ones(1,length(h_bar)))./Cg;
Hs_linear  = sqrt(nansum(16*Snn_linear*df,1));
iLinear = find(abs(imag(Hs_linear))>1e-4,1,'last');
Hs_linear = real(Hs_linear);
Hs_linear(1:iLinear)=nan;
%
% load PVlab pressure sensor data from the 50s cross-shore line (near ring of doom)
load('/home/derek/projects/ShortCrests/ROD/mat_data/RODSEX_09291200_50s_line_pressure.mat','out','rod')
load('/home/derek/projects/ShortCrests/ROD/mat_data/RODSEX_09291200_50s_line_Hs.mat','X','Y','Hs','rodX','rodY','rodHs','ID','P','rodP')
%
% compare crosshore Hs to pressure guages
%
% plot box stuff
xm = 2;
ym = 2;
pw = 8;
ph = 5;
ppos1 = [xm ym pw ph];
ps = [1.2*xm+pw 2*ym+ph];
%
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'color','w','inverthardcopy','off')
%
a1 = axes('units','centimeters','position',ppos1);
p1 = plot(x,Hs_x,'k',x, Hs_linear,':r');hold on
iparos = find(ID=='q');
itrito = find(ID=='p');
p2 = plot(X(itrito),Hs(itrito),'xb',X(iparos),Hs(iparos),'+c',rodX,rodHs,'ob');
set(a1,'xlim',[100 350],'ylim',[0 2],'ytick',[0 0.5 1 1.5 2],'xtick',[100 150 200 250 300 350])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$H_s$ [m]','interpreter','latex')
hh = legend([p1; p2],'modeled','linear shoal','triton','paros','ring-o-doom');
set(hh,'interpreter','latex','fontsize',9,'location','southeast')
fout = sprintf('/home/derek/projects/ShortCrests/MOD/figures/wave_height_0929_%s_v4.png',runID);
exportgraphics(fig,fout)
%
%
% compare mean sealevel elevation
%
h_bar_adv = interp1(x,h_bar,X);
h_bar_rod = interp1(x,h_bar,rodX);
eta_bar_adv = interp1(x,eta_bar,X);
eta_bar_rod = interp1(x,eta_bar,rodX);
%
%
meanSeaLevelDiff = mean([P+h_bar_adv-eta_bar_adv,rodP+h_bar_rod-eta_bar_rod])
sdevSeaLevelDiff = std( [P+h_bar_adv-eta_bar_adv,rodP+h_bar_rod-eta_bar_rod])
%
fig0 = figure('units','centimeters');
pos = get(fig0,'position');
pos(3:4)=ps;
set(fig0,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'color','w','inverthardcopy','off')
%
%
a0 = axes('units','centimeters','position',ppos1);
p1 = plot(x(x>=100),eta_bar(x>=100),'k');hold on
plot([100 x(end)],[0 0],':k')
iparos = find(ID=='q');
itrito = find(ID=='p');
p2 = plot(X(itrito),P(itrito)+h_bar_adv(itrito),'xb',...
          X(iparos),P(iparos)+h_bar_adv(iparos),'+c',rodX,rodP+h_bar_rod,'ob');
set(a0,'xlim',[100 250],'ylim',[-0.5 2],'ytick',[-0.5:0.5:2],'xtick',[100 150 200 250 300 350])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$(\eta,\bar{p}-h)$ [m]','interpreter','latex')
hh = legend([p1; p2],'modeled','triton','paros','ring-o-doom');
set(hh,'interpreter','latex','fontsize',9,'location','southeast')
fout = sprintf('/home/derek/projects/ShortCrests/MOD/figures/mean_sealevel_0929_%s_v4.png',runID);
exportgraphics(fig0,fout)
%
%
% loop over pressure guages and plot Snn at each to compare
I = find(freqs>0.05 & freqs<0.3);
for ii=51:55
    x0 = mean([out(ii).xT,out(ii).xP]);
    if isempty(x0) | isnan(x0), continue, end
    fr = out(ii).fr;
    Ir = find(fr>0.05 & fr<0.3);
    SeT= out(ii).SeT/1e4;
    SeP= out(ii).SeP/1e4;
    Hobs = max([max(SeT),max(SeP)]);% mean([out(ii).HsT, out(ii).HsP]);
    %
    % nearest modeled gridpoint
    [~,ig] = min(abs(x-x0));
    Smod = Snn_x(:,ig);
    Hmod = max(Smod);% 4*sqrt(sum(Smod(I).*df));
    %
    fig = figure('units','centimeters');
    pos = get(fig,'position');
    pos(3:4)=ps;
    set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'color','w','inverthardcopy','off')
    %
    a1 = axes('units','centimeters','position',ppos1);
% $$$     loglog(freqs,Smod,'-k','linewidth',1.5),hold on
    plot(freqs,Smod,'-k','linewidth',1.5),hold on
    lab = {'modeled'};
    if ~isempty(SeT)
% $$$         loglog(fr,SeT,'-b','linewidth',1.5)
        plot(fr,SeT,'-b','linewidth',1.5)
        lab = [lab;{'triton'}];
    end
    if ~isempty(SeP)
% $$$         loglog(fr,SeP,'--b','linewidth',1.5)
        plot(fr,SeP,'--b','linewidth',1.5)
        lab = [lab;{'paros'}];
    end
    title(sprintf('$$ x= %3.0f m$$',x0),'interpreter','latex')
    xlabel('$f$ [Hz]','interpreter','latex')
    ylabel('$S_{\eta\eta}$ [m$^2$/Hz]','interpreter','latex')
% $$$     set(a1,'xtick',[1e-2 1e-1 1e-0],'ytick',[1e-4 1e-3 1e-2 1e-1 1e0],...
% $$$             'yticklabel',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'},'ticklabelinterpreter','latex',...
% $$$             'xlim',[1e-2 1e0],'ylim',[0.9e-4 5e0],'ticklength',2*get(a1,'ticklength'))
    set(a1,'xtick',[0.0 0.1 0.2 0.3],'xminortick','on',...'ytick',[0 0.5 1 1.5 2 2.5],...
            'ticklabelinterpreter','latex','xlim',[0 0.3],'ylim',[0 ceil(max(Hobs,Hmod)*10)/10],'ticklength',2*get(a1,'ticklength'))           
%           'yticklabel',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'},

    a1.XAxis.MinorTickValues = [0.05:0.1:0.3];
    hh = legend(lab);
% $$$     set(hh,'location','southwest','fontsize',9,'interpreter','latex')
    set(hh,'fontsize',9,'interpreter','latex')    
    fout = sprintf('/home/derek/projects/ShortCrests/MOD/figures/spectra_0929_%s_p%2d_v4.png',runID,out(ii).ID);
    exportgraphics(fig,fout)
    %    close(fig)
end
%
%
%
% nearest modeled gridpoint
Irod = find(rod.fr>0.05 & rod.fr<0.3);
Hobs = max(rod.Se(Irod));
[~,ig] = min(abs(x-rod(1).x));
Smod = Snn_x(:,ig);
Hmod = max(Smod);
%
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'color','w','inverthardcopy','off')
%
a1 = axes('units','centimeters','position',ppos1);
% $$$ loglog(freqs,Smod,'-k','linewidth',1.5),hold on
plot(freqs,Smod,'-k','linewidth',1.5),hold on
lab = {'modeled','ring-o-doom'};
plot(rod(1).fr,rod(1).Se,'--b','linewidth',1.5)
% $$$ loglog(rod(1).fr,rod(1).Se,'--b','linewidth',1.5)
title(sprintf('$$ x_\\mathrm{rod}= %3.0f m$$',rod(1).x),'interpreter','latex')
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$S_{\eta\eta}$ [m$^2$/Hz]','interpreter','latex')
% $$$ set(a1,'xtick',[1e-2 1e-1 1e-0],'ytick',[1e-4 1e-3 1e-2 1e-1 1e0],...
% $$$        'yticklabel',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'},'ticklabelinterpreter','latex',...
% $$$        'xlim',[1e-2 1e0],'ylim',[0.9e-4 5e0],'ticklength',2*get(a1,'ticklength'))
set(a1,'xtick',[0.0 0.1 0.2 0.3],'xminortick','on',...'ytick',[0 0.5 1 1.5 2 2.5],...
       'ticklabelinterpreter','latex','xlim',[0 0.3],'ylim',[0 ceil(max(Hobs,Hmod)*10)/10],'ticklength',2*get(a1,'ticklength'))           

hh = legend(lab);
% $$$ set(hh,'location','southwest','fontsize',9,'interpreter','latex')
set(hh,'fontsize',9,'interpreter','latex')    
fout = '/home/derek/projects/ShortCrests/MOD/figures/spectra_0929_run00_ROD_v4.png';
exportgraphics(fig,fout)
% close(fig)
%
%
%
focusTime =   '09/29/2013 12:30:00';
tFocus    = datenum(focusTime);
[t,awacHs,awacTm,awacDm,awacmSpread,awacpSpread,awacfSpread,awacEt,fr,dirs] =load_11mAWAC_stats(2013);
[dtF,itF] = min(abs(tFocus-t)*24*3600);
%
E = awacEt(:,:,itF);%mean(Et(:,:,itF+[-1 0]),3);
awacEf= sum(E.*5,1);
awacHsF = awacHs(itF);
awacTpF= 1./freqs(find(awacEf==max(awacEf)));
awacTmF = awacTm(itF);
om  = freqs*2*pi;
%
[S,~,~]=Jonswap('Omega',2*pi*freqs,'Hs',Hs0,'Tm' ,Tp0, 'Type', 3, 'TEnd',90,'Cap',3,'Gamma',Gam);
% $$$ [S2,~,~]=Jonswap('Omega',2*pi*freqs,'Hs',Hs0,'Tm' ,Tp0, 'Type', 3, 'TEnd',90,'Cap',3,'Gamma',1.5);
%
fig = figure('units','centimeters');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'color','w','inverthardcopy','off')
%
a1 = axes('units','centimeters','position',ppos1);
%
semilogx(freqs,Snn_wg,'-k',fr,awacEf,'-r',freqs,2*pi*S,'--r',rod(1).fr,rod(1).Se,'-b')
set(a1,'xlim',[1/20 1/3])
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$S_{\eta\eta}$ [m$^2$/Hz]','interpreter','latex')
hh = legend('modeled','11m awac','P/M fit','rod');
set(hh,'fontsize',9,'interpreter','latex')
fout = sprintf('/home/derek/projects/ShortCrests/MOD/figures/spectra_0929_%s_offshore_and_ROD_v4.png',runID);
exportgraphics(fig,fout)
% close(fig)

I = find(freqs>0.05 & freqs<0.3);
Hs_mod = 4*sqrt(sum(Snn_wg(I).*df))
I11m = find(fr>0.05 & fr<0.3);
Hs_11m = 4*sqrt(sum(awacEf(I11m).*(fr(2)-fr(1))))
Hs_JONSWAP = 4*sqrt(sum(2*pi*S(I).*df))
Hs_ROD = 4*sqrt(sum(rod.Se(Irod).*(rod.fr(2)-rod.fr(1))))

Tp_mod = 1/freqs(find(Snn_wg==max(Snn_wg)))
Tp_11m = 1/fr(find(S==max(S)))
TP_JONSWAP = 1/freqs(find(awacEf==max(awacEf)))
frod = rod.fr(Irod);
rodE = rod.Se(Irod);
Tp_ROD = 1/frod(find(rodE==max(rodE)))