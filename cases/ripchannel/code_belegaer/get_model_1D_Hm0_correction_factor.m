clear all
close all
%
%
runBATHY = 'planar1D';
[HEIGHT,PERIOD,DIR,SPREAD] = ripchannel_parameter_space(runBATHY);
Hs = [];
for HH=1:length(HEIGHT)
    H = HEIGHT(HH);
    for TT=1:length(PERIOD)
        T = PERIOD(TT);
        for SS=1:length(SPREAD)
            S = SPREAD(SS);               
            for DD=1:length(DIR)
                D = DIR(DD);
                Hs = cat(2,Hs,H);
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
fin = '/data2/ripchannel/mat_data/planar1D_Hm0_error.mat';
load(fin);
%
figure,
x = 1:size(Hm0,2);
plot(x,Hm0'./Hs)
yl = yline(sqrt(2),'--k');
lg = legend(yl,'$\sqrt{2}$')
set(lg,'interpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$H_\mathrm{mod}/H_\mathrm{inp}$ [~]','interpreter','latex')
figName = [info.rootDAT,filesep,runBATHY,'_Hm0_correction.png'];
exportgraphics(gcf,figName)
