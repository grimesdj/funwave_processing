function [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition(u,v,dx,dy);
%function [psi,u_psi,v_psi,phi,u_phi,v_phi]=get_vel_decomposition(u,v,dx,dy);
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

% Make some matrices
[ny,nx]= size(u);
Ly     = (ny-1)*dy;
Lx     = (nx-1)*dx;

% Take some derivitives
ux       = get_dx(u,dx);
ux       = [u(:,1)/dx ux];
uy       = get_dy(u,dy);
vx       = get_dx(v,dx);
vy       = get_dy([v(end,1:end-1);v(1:end-1,1:end-1)],dy);
vbar     = mean(mean(v));
ubar     = mean(mean(u));
psi_at_lx= vbar*Lx;
phi_at_lx= ubar*Lx;

divu         = ux + vy; %forcing for the potential
divu(:,end)  = divu(:,end)-phi_at_lx/dx^2;
curlu        = vx-uy; %forcing for the streamfunction
curlu(:,end) = curlu(:,end)-psi_at_lx/dx^2;

%figure(1)
%subplot(1,2,2)
%pcolor(divu); shading flat
%colorbar('horiz')

% Solves Laplace's equation using FFT in y
un_diag          = 1/dx^2;
alpha            = un_diag*ones(1,nx-1);
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

psi0    = ifft(ifftshift(Xpsi.',1));
u_psi   = -get_dy([psi0(end,:);psi0(1:end-1,:)],dy);
u_psi   = [u_psi(end,:);u_psi(1:end-1,:)];
psi     = [zeros(ny,1) psi0 psi_at_lx*ones(ny,1)];
v_psi   = get_dx(psi,dx);

phi     = ifft(ifftshift(Xphi.',1));
phi     = [phi phi_at_lx*ones(ny,1)];
u_phi   = get_dx(phi,dx);
v_phi   = get_dy(phi,dy);


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
