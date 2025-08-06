function in = time_average_field(in,flt);
%
% USAGE: out = time_average_field(in,flt);
%

sz  = size(in);
in  = reshape( permute(in, [3 1 2]), [sz(3), prod(sz(1:2))] );
in  = conv2(in,flt(:),'same');
in  = permute( reshape( in, [sz(3) sz([1:2])]), [2 3 1]);
