
info.rootSim = '/scratch/grimesdj/ripchannel/flat/output2/'

%% Read final Hs field (Mglob and Nglob from input file)
Mglob = 1364;
Nglob = 2999;
HsFiles = dir([info.rootSim,filesep,'output',filesep,'Hsig_*']);
info.Hs = 0;
Nf      = length(HsFiles);
tmp = 0;
for jj=1:Nf
    fin = [HsFiles(jj).folder,filesep,HsFiles(jj).name];
    fid = fopen(fin);
    dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
    dum = (dum/4).^2;
    tmp = tmp+dum';% so I'm transposing so rows are y-coord and columns are x-coord
    fclose(fid);
end
Hs = 4*sqrt(tmp/Nf);

fid = fopen('/scratch/grimesdj/ripchannel/flat/Hs_OldSponge','w');
fwrite(fid,Hs','*double');% this looks like it's transposed relative to ASCII convension
fclose(fid)


