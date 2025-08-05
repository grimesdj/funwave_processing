function k = wavenumber_FunwaveTVD(om,h);
%
% USAGE: k = wavenumber_FunwaveTVD(om,h);
%
% function that takes the radian wave frequency and
% a vector of depths and returns the wavenumber at that
% depth by solving the linearized dispersion relationship
% for Nwogu 1993's Boussinesq approximation:
%
%     w^2 (h/g) = (kh)^2 (1-a1(kh)^2)/(1-a0(kh)^2);
%
% 
% if 2D, frequency runs down (om needs to be column) and cross-shore depth runs
% across (h needs to be a row)

som=size(om);
sh =size(h);

om=repmat(om,[1 sh(2)]);
h =repmat(h ,[som(1) 1]);

g = 9.81;


a  = -0.39;
a1 = (a + 1/3);
b0 = om.^2.*h/g;
b1 = (1+a.*b0);

k = sqrt(  (b1 - sqrt(b1.*b1-4*a1.*b0) )./(2*a1) )./h;

