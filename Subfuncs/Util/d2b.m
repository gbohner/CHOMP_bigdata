function y = d2b(x,varargin)

%Originally: Copyright (c) 2010, Zacharias Voulgaris
%All rights reserved.
%Fixed by Gergo Bohner 07/03/2016

% Convert a decimanl number into a binary array
% 
% Similar to dec2bin but yields a numerical array instead of a string and is found to
% be rather faster

c = max(ceil(log(x)/log(2)),1); % Number of divisions necessary ( rounding up the log2(x) )
if nargin > 1
  y(max(c,varargin{1})) = 0; % Initialize output array
else
  y(c) = 0;
end
for i = 1:c
    r = floor(x / 2);
    y(c+1-i) = x - 2*r;
    x = r;
end
