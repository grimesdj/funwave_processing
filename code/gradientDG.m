function varargout = gradientDG(f,varargin)
%GRADIENT Approximate gradient.
%
%   This is same as matlab's default except that it does not flip the first and second
%   dimensions
%
%   [FX,FY] = GRADIENT(F) returns the numerical gradient of the matrix F.
%   FX corresponds to dF/dx, the differences in x (horizontal) direction.
%   FY corresponds to dF/dy, the differences in y (vertical) direction.
%   The spacing between points in each direction is assumed to be one.
%   When F is a vector, DF = GRADIENT(F) is the 1-D gradient.
%
%   [FX,FY] = GRADIENT(F,H), where H is a scalar, uses H as the constant
%   spacing between points in each direction.
%
%   [FX,FY] = GRADIENT(F,HX,HY), when F is 2-D, uses the spacing specified
%   by HX and HY. HX and HY can either be scalars to specify the spacing
%   between coordinates or vectors to specify the coordinates of the
%   points. If HX and HY are vectors, their lengths must match the
%   corresponding dimension of F.
%
%   [FX,FY,FZ] = GRADIENT(F), when F is a 3-D array, returns the numerical
%   gradient of F. FZ corresponds to dF/dz, the differences in the z
%   direction. GRADIENT(F,H), where H is a scalar, uses H as the constant
%   spacing between points in each direction.
%
%   [FX,FY,FZ] = GRADIENT(F,HX,HY,HZ) uses the spacing given by HX, HY, HZ.
%
%   [FX,FY,FZ,...] = GRADIENT(F,...) extends similarly when F is N-D and
%   must be invoked with N outputs and either 2 or N+1 inputs.
%
%   Note: The first output FX is always the gradient along the first
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the second dimension of F, going across rows.  For
%   the third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%
%   Examples:
%       x = -2:.2:2; y = (-1:.15:1).';
%       z = x .* exp(-x.^2 - y.^2);
%       [px,py] = gradient(z,.2,.2);
%       contour(x,y,z), hold on
%       quiver(x,y,px,py), hold off
%
%   Class support for input F:
%      float: double, single
%
%   See also DIFF, DEL2.

%   Copyright 1984-2019 The MathWorks, Inc.

[f,ndim,loc,rflag] = parse_inputs(f,varargin);
nargoutchk(0,ndim);

% Loop over each dimension. 
varargout = cell(1,ndim);
siz = size(f);

% first dimension 
prototype = real(full(f([])));
g = zeros(siz, 'like', prototype); % case of singleton dimension
h = loc{1}(:);
n = siz(1);
% Take forward differences on left and right edges
if n > 1
   g(1,:) = (f(2,:) - f(1,:))/(h(2)-h(1));
   g(n,:) = (f(n,:) - f(n-1,:))/(h(end)-h(end-1));
end

% Take centered differences on interior points
if n > 2
   g(2:n-1,:) = (f(3:n,:)-f(1:n-2,:)) ./ (h(3:n) - h(1:n-2));
end

varargout{1} = g;

% second dimensions and beyond
if ndim == 2
    % special case 2-D matrices to support sparse matrices,
    % which lack support for N-D operations including reshape
    % and indexing
    n = siz(2);
    h = reshape(loc{2},1,[]);
    g = zeros(siz, 'like', prototype);
    
    % Take forward differences on left and right edges
    if n > 1
        g(:,1) = (f(:,2) - f(:,1))/(h(2)-h(1));
        g(:,n) = (f(:,n) - f(:,n-1))/(h(end)-h(end-1));
    end
    
    % Take centered differences on interior points
    if n > 2
        h = h(3:n) - h(1:n-2);
        g(:,2:n-1) = (f(:,3:n) - f(:,1:n-2)) ./ h;
    end
    varargout{2} = g;
    
elseif ndim > 2
    % N-D case
    for k = 2:ndim
        n = siz(k);
        newsiz = [prod(siz(1:k-1)) siz(k) prod(siz(k+1:end))];
        nf = reshape(f,newsiz);
        h = reshape(loc{k},1,[]);
        g = zeros(size(nf), 'like', prototype);
        
        % Take forward differences on left and right edges
        if n > 1
            g(:,1,:) = (nf(:,2,:) - nf(:,1,:))/(h(2)-h(1));
            g(:,n,:) = (nf(:,n,:) - nf(:,n-1,:))/(h(end)-h(end-1));
        end
        
        % Take centered differences on interior points
        if n > 2
            h = h(3:n) - h(1:n-2);
            g(:,2:n-1,:) = (nf(:,3:n,:) - nf(:,1:n-2,:)) ./ h;
        end
        
        varargout{k} = reshape(g,siz);
    end
end

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim > 1
    % disp('not swapping 1st and 2nd dims')
    % varargout([2 1]) = varargout([1 2]);
elseif rflag
    varargout{1} = varargout{1}.';
end


%-------------------------------------------------------
function [f,ndim,loc,rowflag] = parse_inputs(f,h)

loc = {}; % spacing along the x,y,z,... directions
ndimsf = ndims(f);
ndim = ndimsf;
rowflag = false;
if isvector(f)
    ndim = 1;
    if isrow(f) % Treat row vector as a column vector
        rowflag = true;
        f = f.';
    end
end

% Default step sizes: hx = hy = hz = 1
indx = size(f);
if isempty(h)
    % gradient(f)
    loc = cell(1,ndimsf);
    for k = 1:ndimsf
        loc(k) = {1:indx(k)};
    end
elseif isscalar(h) % gradient(f,h)
    if isscalar(h{1})
        % Expand scalar step size
        loc = cell(1,ndimsf);
        for k = 1:ndimsf
            loc(k) = {h{1}*(1:indx(k))};
        end
    elseif ndim == 1
        % Check for vector case
        if numel(h{1}) ~= numel(f)
            error(message('MATLAB:gradient:InvalidGridSpacing'));
        end
        loc(1) = h(1);
    else
        error(message('MATLAB:gradient:InvalidInputs'));
    end
elseif ndimsf == numel(h)  % gradient(f,hx,hy,hz,...)
    % Swap 1 and 2 since x is the second dimension and y is the first.
    loc = h;
    if ndim > 1
        % disp('skipping the matlab permute function')
        % loc([2 1]) = loc([1 2]);
    end
    % replace any scalar step-size with corresponding position vector, and
    % check that the values specified in each position vector is the right
    % shape and size
    for k = 1:ndimsf
        if isscalar(loc{k})
            loc{k} = loc{k}*(1:indx(k));
        elseif ~isvector(squeeze(loc{k})) || numel(loc{k}) ~= size(f, k)
            error(message('MATLAB:gradient:InvalidGridSpacing'));
        end
    end 
else
    error(message('MATLAB:gradient:InvalidInputs'));
end
