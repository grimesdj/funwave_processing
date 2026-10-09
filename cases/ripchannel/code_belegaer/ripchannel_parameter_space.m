function [HEIGHT,PERIOD,DIR,SPREAD,runWAVESlist,grids] = ripchannel_parameter_space(runBATHY);
%
% USAGE: [HEIGHT,PERIOD,DIR,SPREAD,runWAVESlist,GRIDS] = ripchannel_parameter_space(runBATHY);
%
% e.g., runBATHY='planar1D' outputs DIR=0, SPREAD=0

grids = {};
switch runBATHY
    %% SUMMER 2026 GRIDS/RUNS
  case {'gapRip'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 10];
    grids  = {'barRip2','barRip3','terRip0','terRip2','terRip3'};
  case {'gapRip1'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 10];
    grids  = {'bar1Rip0','bar1Rip1','bar1Rip2','bar1Rip3','ter1Rip0','ter1Rip1','ter1Rip2','ter1Rip3'};
  case {'spreadRip1'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 2 4 10 20];
    grids  = {'bar1Rip1','ter1Rip1'};
  case {'highRipXtraGaps'}
    HEIGHT = [1.5];
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [10];
    grids  = {'barRip2','barRip3','terRip0','terRip2','terRip3'};
  case {'highRip1'}
    HEIGHT = [0.5 1.0 1.5];
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 10];
    grids  = {'bar1Rip1','ter1Rip1'};
  case {'highRip1XtraGaps'}
    HEIGHT = [0.5 1.5];
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [10];
    grids  = {'bar1Rip0','bar1Rip2','bar1Rip3','ter1Rip0','ter1Rip2','ter1Rip3'};
  %% OSM2026 GRIDS/RUNS: NB, the 0-degree spread runs are being converted to 1-deg; and highRip 0.5m run is removed.  
  case {'testRip'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = 10;
    grids  = {'barRip0','barRip1','terRip1'};
  case {'spreadRip'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 2 4 10 20];
    grids  = {'barRip1','terRip1'};
  case {'uniRip'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [2 4 10 20];
    grids  = {'bar2D','ter2D'};
  case {'highRip'}
    HEIGHT = [1.5];
    PERIOD = 10;
    DIR    = 0;
    SPREAD = [1 10];
    grids  = {'barRip1','terRip1'};
  %% TESTING GRIDS/RUNS
  case {'planar2D','barred2D','terraced2D'}
    HEIGHT = [0.5 1 1.5];
    PERIOD = [8 10];
    DIR    = 0;
    SPREAD = [0 5 10 20];
  case {'test2D'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = 10;
    grids  = {'plnr2D','bar2D','ter2D'};
  case {'planar1D','barred1D','terraced1D','barred1Ddx050','barred1Ddx025'}
    HEIGHT = [0.5 1 1.5];
    PERIOD = [8 10];
    DIR    = 0;
    SPREAD = 0;
  case {'resolution1D'}
    HEIGHT = 1;
    PERIOD = 10;
    DIR    = 0;
    SPREAD = 0;
    grids  = {'terraced1Ddx025','terraced1Ddx050','terraced1D','planar1Ddx025','planar1Ddx050','planar1D'};
end


idx=1;
for HH=1:length(HEIGHT)
    H = HEIGHT(HH);
    for TT=1:length(PERIOD)
        T = PERIOD(TT);
        for SS=1:length(SPREAD)
            S = SPREAD(SS);               
            for DD=1:length(DIR)
                D = DIR(DD);
                %
                % construct identifier
                runWAVES = sprintf('h%02dt%02ds%02dd%02d',H*10,T,S,D);
                runWAVESlist{idx} = runWAVES;
                idx=idx+1;
            end
        end
    end
end
