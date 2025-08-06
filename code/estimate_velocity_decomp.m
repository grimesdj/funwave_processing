function [Urot,Vrot,VORT] = estimate_velocity_decomp(u,v)
%
%
%

[omega,~] = curl(xx,yy,U,V);
