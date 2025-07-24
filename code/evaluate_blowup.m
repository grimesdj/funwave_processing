% last waterlevel file (ASCII)
fin = '/home/derek/projects/ShortCrests/MOD/funwave/0929/run03/output/eta_99999';
% grid specifications
% Mglob=Nx-1;
% Nglob=Ny-1;
Mglob = 1700;
Nglob = 1500;
dx    = 0.5;
dy    = 1.0;
% peak period of waves
Tp    = 10;
% wave maker
xWM   = 692;
hWM   = 7;
dWM   = 0.5;
% sponge regions
sew   = 80;
sww   = 5;
%
%
x = 0.5*[0:Mglob-1];
y = 1*[0:Nglob-1];
%
eta = load(fin);
figure, imagesc(x,y,eta),colormap(cmocean('ice')),caxis([-1.5 1.5])

xline(xWM,'-r','linewidth',2)

% estimate limits on WM region
k = wavenumber(2*pi/Tp,hWM);
L = 2*pi/k;
W = dWM*L/2;

xline(xWM + W/2*[-1 1],'--r','linewidth',2)
xline([sww x(end)-sew],'--g','linewidth',2)