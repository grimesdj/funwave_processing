function [Fbx,Fby] = estimate_Fbr(p,q,nubrk,H,dx,dy)
%
% USAGE: [Fbx,Fby] = estimate_Fbr(p,q,nubrk,H,dx,dy);
%
% estimate breaking wave forcing

% water depth
% $$$ H        = h+eta;
H(H<=0)  = 0.001;
% gradient of fluxes
[Py, Px]   = gradientDG(p);% grad(H.*u);
[Qy, Qx]   = gradientDG(q);% grad(H.*v);
% divergence of fluxes... stress
[Pyy,   ~] = gradientDG(nubrk.*Py/dy);
[~  , Pxx] = gradientDG(nubrk.*Px/dx);
[Qyy,   ~] = gradientDG(nubrk.*Qy/dy);
[~  , Qxx] = gradientDG(nubrk.*Qx/dx);
%
clear Py Px Qy Qx
Pyy           = Pyy/dy; 
Pxx           = Pxx/dx;
Qyy           = Qyy/dy;
Qxx           = Qxx/dx;
% body forcing / per unit volume
Fbx             = (Pyy+Pxx)./H; clear Pxx Pyy
Fby             = (Qxx+Qyy)./H; clear Qxx Qyy
 
