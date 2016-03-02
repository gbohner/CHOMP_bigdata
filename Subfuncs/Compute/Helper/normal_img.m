function [y, A, B] = normal_img(I, sig1, sig2, varargin)
% I = mean image, sig 1 and 2 are the filter sizes for mean and contrast
% normalization, varargin{1} can be the temporal variance, which is used in
% stead of spatial variance then.
y = I;

if sig1 == 0
  A = ones(size(y));
  B = ones(size(y));
  return;
end

Mask = ones(size(I));

lx = ceil(2*sig2);
dt = (-lx:lx)';

% keyboard;
if sig1<.1 %Hacky solution to get rid of line artefact
  filter = zeros(length(dt));
  filter(lx+1, :) = 1;
elseif sig1<.25
    filter = zeros(length(dt));
    filter(lx:lx+2, lx:lx+2) = 1;
else
    sig = sig1;
    filter = exp(-dt.^2/(2*sig^2)) * exp(-dt'.^2/(2*sig.^2));
end

filter = filter/sum(filter(:));
Norms = conv2(Mask, filter, 'same');

A = conv2(y, filter, 'same');
A = A./Norms;

y = (y - A) ;

sig = sig2;
filter = exp(-dt.^2/(2*sig^2)) * exp(-dt'.^2/(2*sig.^2));
filter = filter/sum(filter(:));
Norms = conv2(Mask, filter, 'same');

if nargin == 3
  B= conv2(y.^2, filter, 'same');
elseif nargin == 4
  B= conv2(varargin{1}, filter, 'same'); %smooth the true variance
end
B = B./Norms;

y = y./B.^.5;



