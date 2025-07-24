function info = advance_FUNWAVE_wave_input_phase(info,input_files,final_times,output_identifiers);
%
%
%
for jj = 1:length(input_files)
    fin = input_files{jj};
    t   = final_times(jj);
    %
    %
    fid = fopen(fin);
    str = fgetl(fid);
    str = split(str,'-');
    Nfrq = str2num(str{1});
    %
    str = fgetl(fid);
    str = split(str,'-');
    Tp  = str2num(str{1});
    %
    for kk = 1:Nfrq
        str  = fgetl(fid);
        str  = split(str,'- Freq');
        f(kk,1)= str2num(str{1});
    end
    %
    for kk = 1:Nfrq
        str  = fgetl(fid);
        str  = split(str,'- Dire');
        d(kk,1)= str2num(str{1});
    end
    %
    for kk = 1:Nfrq
        str  = fgetl(fid);
        a(kk,1)= str2num(str);
    end
    %
    for kk = 1:Nfrq
        str  = fgetl(fid);
        p(kk,1)= str2num(str);
    end
    %
    % advance the phase
    p_new = p+360*f*t;
    p_new = mod(p_new,360);
    %
    fileStr = split(fin,'.');
    fout    = sprintf([fileStr{1},'_',output_identifiers{jj},'_advanced_%ds.',fileStr{2}],round(t));
    %
    fid       = fopen(fout,'w');
    fprintf(fid,'%5i      - NumFreq  \n',Nfrq);
    fprintf(fid,'%10.3f   - PeakPeriod \n',Tp);
    fprintf(fid,'%2.8f   - Freq \n',f);
    fprintf(fid,'%2.8f   - Dire \n',d);
    fclose(fid);
    dlmwrite(fout,a,'delimiter','\t','-append','precision',5);
    dlmwrite(fout,p_new,'delimiter','\t','-append','precision',5);
end