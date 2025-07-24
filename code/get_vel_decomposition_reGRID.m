function [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition_reGRID(u,v,dx,dy);
% function [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition(u,v,dx,dy);
% This function returns the velocity stream function and velocity potential
% given a velocity field (u,v). We assume zero velocity on the x
% boundaries. The equation being solved here is fo the form:
% ui^+vj^=div(phi)+curl(psi)
% INPUT:
% u   = cross-shore velocity (ms^{-1})
% v   = alongshore velocity (ms^{-1})
% dx  = cross-shore grid spacing (m)
% dy  = alongshore grid spacing (m)
% OUTPUT:
% psi = Velocity stream function
% phi = Velocity potential
%
% This function has been written by Dr. Matthew Spydell, SIO, UCSD
% edited for regular grid (non-C-grid) by Derek Grimes, UNCW

% Make some matrices
[ny,nx]= size(u);
% Ly     = (ny-1)*dy;
% Lx     = (nx-1)*dx;
Ly = ny*dy;
Lx = nx*dx;

% Convert nan-mask to zeros
inan    = isnan(u) | isnan(v);
u(inan) = 0;
v(inan) = 0;

% $$$ % enforce (u,v)=0 on x-boundaries? 
% $$$ u(:,[1 nx]) = 0;
% $$$ v(:,[1 nx]) = 0;

% Take some derivitives, gradientDG does not permute rows/columns 1-2
[uy, ux]  = gradientDG(u,dy,dx);
[vy, vx]  = gradientDG(v,dy,dx);
vbar     = mean(mean(v));
ubar     = mean(mean(u));
psi_at_lx= vbar*Lx;
phi_at_lx= ubar*Lx;

divu         = ux + vy; %forcing for the potential
divu(:,nx)  = divu(:,nx)-phi_at_lx/dx^2;
curlu        = vx-uy; %forcing for the streamfunction
curlu(:,nx) = curlu(:,nx)-psi_at_lx/dx^2;

%figure(1)
%subplot(1,2,2)
%pcolor(divu); shading flat
%colorbar('horiz')

% Solves Laplace's equation using FFT in y
un_diag          = 1/dx^2;
% alpha            = un_diag*ones(1,nx-1);
alpha            = un_diag*ones(1,nx);
on_diag          = -2/dx^2;
beta             = on_diag*ones(1,nx);
beta_noflux      = beta;
beta_noflux(end) = 2/dx^2;
ov_diag          = 1/dx^2;
gamma            = alpha;
gamma_noflux     = gamma;
gamma_noflux(1)  = 2/dx^2;

Gpsi             = fft(curlu);
Gpsi             = fftshift(Gpsi,1);
Gphi             = fft(divu);
Gphi             = fftshift(Gphi,1);

N                = [-ny/2:(ny/2)-1].';
Xpsi             = zeros(size(Gpsi.'));
Xphi             = zeros(size(Gphi.'));

for a=1:ny
    c         = -4*pi^2*N(a)^2/Ly^2; %it should be -4... 
    g         = Gpsi(a,:).';
    Xpsi(:,a) = tridiag_mat(alpha,beta+c,gamma,g);
    g         = Gphi(a,:).';
    Xphi(:,a) = tridiag_mat(alpha,beta+c,gamma_noflux,g);
end

psi0        = real( ifft(ifftshift(Xpsi.',1)) );
[u_psi,~]   = gradientDG([psi0(end,:);psi0(1:end-1,:)],dy,dx);
u_psi       =-[u_psi(end,:); u_psi(1:end-1,:)];
psi         = [zeros(ny,1) psi0 psi_at_lx*ones(ny,1)];
[~,v_psi]   = gradientDG(psi,dy,dx);
v_psi = v_psi(:,2:end-1);
psi         = psi(:,2:end-1);

phi         = real( ifft(ifftshift(Xphi.',1)) );
phi         = [phi phi_at_lx*ones(ny,1)];
[~,u_phi]   = gradientDG(phi,dy,dx);
u_phi=u_phi(:,2:end);
[v_phi,~]   = gradientDG(phi,dy,dx);
v_phi=v_phi(:,2:end);
phi         = phi(:,2:end);

function dxf=get_dx(f,dx)
% on eta grid!!!
[ny,nx]=size(f);
dxf=1/dx*(f(:,2:end)-f(:,1:end-1));
%dxf=[.5/dx*f(:,2) dxf -.5/dx*f(:,nx-1)];


function dyf=get_dy(f,dy)
[ny,nx]=size(f);
f1=f(2:end,:);
f1=[f1;f(1,:)];
f2=f;
dyf=1/dy*(f1-f2);

function x2=tridiag_mat(alpha,beta,gamma,b)
% Solve the tridiagonal system Ax=b, alpha is under the diagonal,  beta is the
% diagonal, and gamma is above the diagonal

N=length(b);

% Perform forward elimination
for i=2:N
    coeff = alpha(i-1)/beta(i-1);
    beta(i) = beta(i) - coeff*gamma(i-1);
    b(i) = b(i) - coeff*b(i-1);
end

% Perform back substitution
x2(N) = b(N)/beta(N);
for i=N-1:-1:1
    x2(i) = (b(i) - gamma(i)*x2(i+1))/beta(i);
end
x2 = x2.';
