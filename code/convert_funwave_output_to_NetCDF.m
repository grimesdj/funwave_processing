function flog = convert_funwave_output_to_NetCDF(rin,rout,vars,t0,dt0,dx,spanx,dy,spany,rmfiles,Nt,tlimits,isBINARY,Mglob,Nglob)
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
    Nt=inf;
end
%
if ~exist('tlimits','var')
    tlimits = [-inf inf];
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
    fprintf('removing old files:\n')
    dir([rout,var,'_*.nc'])
    eval(['!rm ',rout,var,'*.nc'])
    %
    if (nf~=1) & (nf~=nt)
        disp(['time vector is off: nt=',num2str(nt),'and nf=',num2str(nf)])
        %
        disp('shortenning records to match')
        nf = min(nf,nt);
        nt = nf;
        t0  = t0(1:nt);
        dt0 = dt0(1:nt);
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
        fin = [rin,files(kk).name];
        if ( (kk==1 | flag) & ~isBINARY )
            % read in first file to get domain dimensions
            dum = load(fin);
            [ny,nx] = size(dum);
            %
            x   = [0:nx-1]*dx;
            y   = [0:ny-1]'*dy;
            % if only partially saving variables...
            nx0 = floor(nx/spanx);
            ny0 = floor(ny/spany);
            % pre-allocate, unless this is a single file variable
            if nf==1
                eval([var,'=(dum(1:spany:ny0*spany,1:spanx:nx0*spanx));']),
            else
                eval([var,'=nan(ny0,nx0,min(nt-Nk*Nt,Nt));'])
                % first element is southwest-side of domain
                eval([var,'(:,:,1)=(dum(1:spany:ny0*spany,1:spanx:nx0*spanx));']),
            end
            x = x(1:spanx:nx0*spanx);
            y = y(1:spany:ny0*spany);
            %
            clear dum
            flag=0;
            continue
        end
        %
        if isBINARY
            fid = fopen(fin);
            dum = fread(fid,[Mglob Nglob],'*double');% this looks like it's transposed relative to ASCII convension
            dum = dum';% so I'm transposing so rows are y-coord and columns are x-coord
            fclose(fid);
        else
            dum = load(fin);
        end
        eval([var,'(:,:,kk-Nk*Nt)=(dum(1:spany:ny,1:spanx:nx));']),
        if (kk-Nk*Nt)==Nt
            flag=1;
            Nk  = Nk+1;
            Nj  = Nj+1;
            % get current time segment
            t = t0((Nk-1)*Nt+1:Nk*Nt);
            dt=dt0((Nk-1)*Nt+1:Nk*Nt);
            % partial save for large files
            fout = sprintf([rout,var,'_%02d.mat'],Nk);
            fprintf(' archiving: %s \n',fout)
            nccreate(fout,var,'Format','netcdf4')
            eval(['ncwrite(fout,''',var,''',var);'])
            nccreate(fout,'x','Format','netcdf4')
            ncwrite (fout,'x',x);
            nccreate(fout,'y','Format','netcdf4')
            ncwrite (fout,'y',y);
            nccreate(fout,'t','Format','netcdf4')
            ncwrite (fout,'t',t);
            nccreate(fout,'dt','Format','netcdf4')
            ncwrite (fout,'dt',dt);
            % save('-v7.3',fout,var,'t','dt')
            eval(['clear ',var])            
            flog{Nj,1} = fout;
        end
    end
    %
    if  isinf(Nt) || Nk==0
        % get current time segment
        t = t0;
        dt=dt0;
        fout = ([rout,var,'.mat']);
    elseif nt-(Nk*Nt)>0.25*Nt
        Nk = Nk+1;
        % get current time segment
        t = t0((Nk-1)*Nt+1:nt);
        dt=dt0((Nk-1)*Nt+1:nt);
        fout = sprintf([rout,var,'_%02d.mat'],Nk);
    else
        eval(['clear ',var])            
    end
    %
    if exist(var,'var')
        fprintf(' archiving: %s \n',fout)
        nccreate(fout,var,'Format','netcdf4')
        eval(['ncwrite(fout,''',var,''',var);'])
        nccreate(fout,'x','Format','netcdf4')
        ncwrite (fout,'x',x);
        nccreate(fout,'y','Format','netcdf4')
        ncwrite (fout,'y',y);
        nccreate(fout,'t','Format','netcdf4')
        ncwrite (fout,'t',t);
        nccreate(fout,'dt','Format','netcdf4')
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
