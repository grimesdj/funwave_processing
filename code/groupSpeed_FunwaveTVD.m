function cg = groupSpeed_FunwaveTVD(om,h)
%
% USAGE: cg = groupSpeed_FunwaveTVD(om,h);
%
% Function that takes the radian wave frequency and a vector of
% depths and returns the group velocity, d(omega)/dk, using the
% linearized dispersion relationship for Nwogu (1993)'s
% Boussinesq approximation:
%
%     w^2 (h/g) = (kh)^2 (1-a1(kh)^2)/(1-a(kh)^2)
%
% where
%
%     a  = -0.39
%     a1 = a + 1/3
%
% If 2D, frequency runs down (om needs to be a column) and
% cross-shore depth runs across (h needs to be a row).
%
% INPUT:
%   om  = radian wave frequency [rad/s]
%   h   = water depth [m]
%
% OUTPUT:
%   cg  = group velocity d(omega)/dk [m/s]
%

% Match dimensions used by wavenumber_FunwaveTVD
som = size(om);
sh  = size(h);

om = repmat(om,[1 sh(2)]);
h  = repmat(h,[som(1) 1]);

g = 9.81;

% Nwogu/FUNWAVE coefficients
a  = -0.39;
a1 = a + 1/3;

% -------------------------------------------------------------
% Calculate wavenumber using the same analytical inversion
% as wavenumber_FunwaveTVD.m
% -------------------------------------------------------------
b0 = om.^2 .* h ./ g;
b1 = 1 + a.*b0;

k = sqrt( (b1 - sqrt(b1.^2 - 4*a1.*b0)) ./ (2*a1) ) ./ h;

% Nondimensional wavenumber
kh = k.*h;

% -------------------------------------------------------------
% Group velocity
%
% omega^2 = (g/h) F(kh)
%
% where
%
% F(x) = x^2 (1-a1*x^2)/(1-a*x^2)
%
% Therefore:
%
% cg = d(omega)/dk
%    = g/(2*omega) * dF/dx
%
% -------------------------------------------------------------

x = kh;

N  = x.^2 - a1.*x.^4;
D  = 1 - a.*x.^2;

dN = 2.*x - 4.*a1.*x.^3;
dD = -2.*a.*x;

dFdx = (dN.*D - N.*dD)./D.^2;

cg = g./(2.*om) .* dFdx;

end
