function out = collect_inp( varargin )
%COLLECT_INP Summary of this function goes here
%   Detailed explanation goes here
  out = zeros(nargin,1);

  for i1 = 1:nargin
    out(i1) = varargin{i1};
  end

end

