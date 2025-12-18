clear all
close all
% extract the wavemaker E(f,d) info from log file:
fin = '~/Grants/NSF_SED/ripchannel/test2Davg/LOG.txt';
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
elseif ismember('Input Dire',tmp2)
    iter=iter+1;
    line=fgetl(fid);
    num = str2num(line);
    freq=[freq,num(1)];
    dire0=[dire0,num(2)];
    dire=[dire,num(3)];
    if length(num)==3
        line=fgetl(fid);
        num = str2num(line);
        amp = [amp,num(1)];
        phi = [phi,num(2)];
    elseif length(num)==5
        amp = [amp,num(4)];
        phi = [phi,num(5)];
    end
end

end