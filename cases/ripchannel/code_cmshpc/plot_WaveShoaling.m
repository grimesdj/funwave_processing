spread=0;
dat.height=1;
dat.direction=0;
dat.period=10;
jj=1
fin = '/data2/ripchannel/spreadRip/mat_data/funwave_barRip1_h10t10s00d00_MomentumTerms.nc';
Hs = ncread(fin,'Hsig');
Hs = mean(Hs,3);
if spread==0
    fltWave = hamming(25); fltWave = fltWave'/sum(fltWave);
    Hs = conv2(Hs,fltWave,'same');
end
%
NAME = 'spreadRip-barRip1';
outDIR   = '/data2/ripchannel/mat_data/'
figDIR   = '/data2/ripchannel/figures/'
dat = load([outDIR,'BulkVelocityStats_',NAME,'.mat']);
        %% Plot Moulton's scaling:
        g = 9.8;
        gam = 0.35;
        %% Kludge to find bar-crest:
        icrest = find(dat.x>=dat.bar_location(jj),1,'first');
        hcrest = 0.03*(dat.x(icrest)-50)-dat.bar_amplitude(jj).*exp(-0.5.*((dat.x(icrest)-dat.bar_location(jj))./dat.bar_width(jj)).^2);
        %
        %% attempt 1 to estimate Hbr:
        k0 = wavenumber_FunwaveTVD(2*pi/dat.period(jj),9);
        disper0 = sqrt(g.*k0.*tanh(k0.*9));
        dwdk_full0 = 0.5*(g*tanh(k0.*9) + g*k0.*9.*sech(k0.*9).^2)./disper0;
        Hbr0 = nthroot( dat.height(jj).^4.*dwdk_full0.^2.*(gam/g), 5);
        %% attempt 2 to estimate Hbr
        k  = wavenumber_FunwaveTVD(2*pi/dat.period(jj),dat.h0);
        disper = sqrt(g.*k.*tanh(k.*dat.h0));
        dwdk_full = 0.5*(g*tanh(k.*dat.h0) + g*k.*dat.h0.*sech(k.*dat.h0).^2)./disper;
        Hshoal = real(sqrt( dat.height(jj).^2.*dwdk_full(end)./dwdk_full ));
        Hmax   = 0.88./k.*tanh(gam*k.*dat.h0/0.88);
        ibrk   = find(Hshoal>=Hmax,1,'last');
        Hbr = Hmax(ibrk);
        dETA_moulton0 = -gam.^2/16.*(cosd(dat.direction(jj)).^2 + 0.5).*max(Hbr0/gam-hcrest,0);
        dETA_moulton1 = -gam.^2/16.*(cosd(dat.direction(jj)).^2 + 0.5).*max(Hbr/gam-hcrest ,0);       
        Uscale0(jj) = sqrt(-2*g*dETA_moulton0);
        Uscale1(jj) = sqrt(-2*g*dETA_moulton1);       


        figure, plot(dat.x,Hs(1,:),'-k',dat.x,Hmax,'--r',dat.x,Hshoal,'--b','linewidth',2)
        hold on,yline(Hbr0,'--g','linewidth',2)
        xlabel('$x$','interpreter','latex')
        legend({'H_s','H_{max}','H_{shoal}','H_{shallow}'},'location','northwest')
        set(gca,'xlim',[50 500], 'ylim',[0 4],'tickdir','out')
        exportgraphics(gcf,[figDIR,filesep,'shoal_Hs_to_breaking.pdf'])