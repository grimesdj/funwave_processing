% clear all
% close all
% extract the wavemaker E(f,d) info from log file:
% fin = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/SED/ripchannel/spreadRip/barRip0_h10t10s02d00/LOG.txt';
% fin = '/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/SED/ripchannel/spreadRip/barRip0_h10t10s00d00/LOG.txt';
function [freq, dire, amp] = parse_wavemaker_info_LOG_aggregated(fin,plotter);
fid = fopen(fin);
iter = 0;

% number of frequencies/directions
Nf_flag = 0;
Nt_flag = 0;

freq = [];
dire0= [];
dire = [];
amp  = [];
phi  = [];
while ~feof(fid)
iter = iter+1;    
line = fgetl(fid);

tmp1 = split(line,'=');
tmp1 = strtrim(tmp1);

tmp2 = split(line,',');

if ismember('Nfreq',deblank(tmp1)) & ~Nf_flag
    Nf = str2num(tmp1{2});
    Nf_flag=1;
elseif ismember('Ntheta',deblank(tmp1)) & ~Nt_flag
    Nt = str2num(tmp1{2});
    Nt_flag=1;
elseif ismember('PBC Dire',tmp2)
    iter=iter+1;
    line=fgetl(fid);
    num = str2num(line);
    freq=[freq,num(1)];
    dire=[dire,num(2)];
    amp =[amp, num(3)];
    phi =[phi, num(4)];
end

end

if plotter
    figure,
    scatter(freq,dire,20,amp,'filled')
end

end

