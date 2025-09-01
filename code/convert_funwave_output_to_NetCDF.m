function flog = convert_funwave_output_to_NetCDF(rin,rout,vars,t0,dt0,dx,spanx,rngx,dy,spany,rngy,rmfiles,Nt,tlimits,isBINARY,Mglob,Nglob)
%
% USAGE: fout = convert_funwave_output_to_NetCDF(root_dir,vars,t,dt,rmfiles)
% need to include t, dt, spanx=span in columns, spany=span in rows; Nt is number of time-steps to break record into.
% neet to make this binary/ascii compatible
nv = length(vars);
nt = length(t0);
%
% rin = rootOut; rout = [rootMat,rootName];t0=t;dt0=dt;
%
if ~exist('Nt','var')
    timeFlag = 1;
    Nt=inf
else
    timeFlag = 0;
end
%
if ~exist('tlimits','var')
    tlimits = [-inf inf];
end
%
rngxFlag=0;
rngyFlag=0;
if any(isinf(rngx))
    rngxFlag=1;
end
if any(isinf(rngy))
    rngyFlag=1;
end
%
if ~exist('isBINARY','var')
    isBINARY=0;
    if isBINARY & ~exist('Mglob','var')
        disp('error using binary format, need Mglob/Nglob')
        return
    end
end
%
% $$$ if ~exist('Ix','var')
% $$$     Ixs = '1' ;
% $$$     Ixe = 'nx';
% $$$ else
% $$$     Ixs = num2str(Ix(1));
% $$$     Ixe = num2str(Ix(2));
% $$$ end
% $$$ %
% $$$ if ~exist('Iy','var')
% $$$     Iys = '1' ;
% $$$     Iye = 'ny';
% $$$ else
% $$$     Iys = num2str(Iy(1));
% $$$     Iye = num2str(Iy(2));
% $$$ end
%
% iterate over vars
Nj = 0;
for jj = 1:nv
    var = vars{jj};
    nl  = length(var);
    %
    % get struct-array w/ filenames
    files = dir([rin,var,'*']);
    nf    = length(files);
    %
    disp(['working on: ', var])
    % remove old files from archive
    fprintf('\tremoving old files:\n')
    if nf==1 | Nt==inf
        eval(['!rm ',rout,var,'*.nc'])
    else
        eval(['!rm ',rout,var,'_*.nc'])
    end
    %
    if (nf~=1) & (nf~=nt)
        disp(['time vector is off: nt=',num2str(nt),'and nf=',num2str(nf)])
        %
        disp('shortenning records to match')
        nf = min(nf,nt);
        nt = nf;
        t0  = t0(1:nt);
        dt0 = dt0(1:nt);
% $$$         if timeFlag;
% $$$             Nt=nt;
% $$$         end
    end
    %
% $$$     % if files aren't organized... need index of time coordinate
% $$$     if nf>1
% $$$         lead0 = length(files(1).name)-nl+1;
% $$$     end
    %
    flag=0;
    Nk  =0;
    for kk = 1:nf
        iter = Nk*Nt; iter(isnan(iter))=0;
        fin = [rin,files(kk).name];
        if ( (kk==1 | flag) )
            if ~isBINARY
            % read in first file to get domain dimensions
                dum = load(fin,'-ascii');
            else
                fid = fopen(fin);
                dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
                dum = dum';% so I'm transposing so rows are y-coord and columns are x-coord
                fclose(fid);
            end
            [ny,nx] = size(dum);
            % construct (x,y) coordinates
            x   = [0:nx-1]*dx;
            y   = [0:ny-1]'*dy;
            %  define limits on (x,y)
            if rngxFlag, rngx = [1 nx]; rngxFlag=0; end
            if rngyFlag, rngy = [1 ny]; rngyFlag=0; end
            % if only partially saving variables...
            nx0 = floor(diff(rngx)/spanx);
            ny0 = floor(diff(rngy)/spany);
            % pre-allocate, unless this is a single file variable
            if nf==1 
                eval([var,'=(dum(rngy(1):spany:ny0*spany,rngx(1):spanx:nx0*spanx));']),
            else
                eval([var,'=nan(ny0,nx0,min(nt-iter,Nt));'])
                % first element is southwest-side of domain
                eval([var,'(:,:,1)=(dum(rngy(1):spany:ny0*spany,rngx(1):spanx:nx0*spanx));']),
            end
            x = x(1:spanx:nx0*spanx);
            y = y(1:spany:ny0*spany);
            %
            clear dum
            flag=0;
            continue
        end
        % need to incorporate this into above if statement
        if isBINARY
            fid = fopen(fin);
            dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
            dum = dum';% so I'm transposing so rows are y-coord and columns are x-coord
            fclose(fid);
        else
            dum = load(fin,'-ascii');
        end
        %
        eval([var,'(:,:,kk-iter)=(dum(rngy(1):spany:ny0*spany,rngx(1):spanx:nx0*spanx));']),
        %
        if (kk-iter)==min(nt,Nt)
            flag=1;
            Nk  = Nk+1;
            Nj  = Nj+1;
            % get current time segment
            iter0 = (Nk-1)*Nt; iter0(isnan(iter0))=0; 
            t = t0(iter0+1:min(nt,Nk*Nt));
            dt=dt0(iter0+1:min(nt,Nk*Nt));
            % partial save for large files
            fout = sprintf([rout,var,'_%02d.nc'],Nk);
            fprintf(' archiving: %s \n',fout)
            if length(size(eval(var)))==2
                dim = {"y",length(y),"x",length(x)};
            elseif length(size(eval(var)))==3
                dim = {"y",length(y),"x",length(x),"t",length(t)};
            end
            nccreate(fout,var,'Dimensions',dim,'Format','netcdf4')
            eval(['ncwrite(fout,''',var,''',',var,');'])
            nccreate(fout,'x','Dimensions',{"x",length(x)},'Format','netcdf4')        
            ncwrite (fout,'x',x);%      
            nccreate(fout,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
            ncwrite (fout,'y',y);
            nccreate(fout,'t','Dimensions',{"t",length(t)},'Format','netcdf4')
            ncwrite (fout,'t',t);
            nccreate(fout,'dt','Dimensions',{"t",length(t)},'Format','netcdf4')
            ncwrite (fout,'dt',dt);
            % save('-v7.3',fout,var,'t','dt')
            eval(['clear ',var])            
            flog{Nj,1} = fout;
        end
    end
    %
    if  Nk==0
        % get current time segment
        t = t0;
        dt=dt0;
        fout = ([rout,var,'.nc']);
    elseif (nt-iter>0.25*Nt) 
        Nk = Nk+1;
        % get current time segment
        t = t0((Nk-1)*Nt+1:nt);
        dt=dt0((Nk-1)*Nt+1:nt);
        fout = sprintf([rout,var,'_%02d.nc'],Nk);
    else
        eval(['clear ',var])            
    end
    %
    if exist(var,'var')
        fprintf(' archiving: %s \n',fout)
        if length(size(eval(var)))==2
            dim = {"y",length(y),"x",length(x)};
        elseif length(size(eval(var)))==3
            dim = {"y",length(y),"x",length(x),"t",length(t)};
        end
        nccreate(fout,var,'Dimensions',dim,'Format','netcdf4')
        eval(['ncwrite(fout,''',var,''',',var,');'])
        nccreate(fout,'x','Dimensions',{"x",length(x)},'Format','netcdf4')        
        ncwrite (fout,'x',x);%      
        nccreate(fout,'y','Dimensions',{"y",length(y)},'Format','netcdf4')
        ncwrite (fout,'y',y);
        nccreate(fout,'t','Dimensions',{"t",length(t)},'Format','netcdf4')
        ncwrite (fout,'t',t);
        nccreate(fout,'dt','Dimensions',{"t",length(t)},'Format','netcdf4')
        ncwrite (fout,'dt',dt);
% $$$         save('-v7.3',fout,var,'t','dt')
        eval(['clear ',var])
    end
    %
    % remove originals?
    if rmfiles
        if nf>1
            eval(['!rm ',rin,var,'_*'])
        else
            eval(['!rm ',rin,var,'.out'])
        end
    end
    %
    Nj = Nj+1;
    flog{Nj,1} = fout;
end
