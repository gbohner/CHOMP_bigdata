function C = mply( A, B, dims )
%MPLY A B are arbitrary tensors, with the last dimension of A and first
%dimension of B being the same
%   With dims you can specify how many dimensions at the end of A do you
%   want to merge with leading dimensions of B
%
%   Written by Gergo Bohner <gbohner@gatsby.ucl.ac.uk> - 2015/07/13

if nargin <3 
  dims = 1;
end

szA = size(A);
szB = size(B);

if dims == 0
  %kind of kronecker product, just replicate A matrix till it matches the
  %dimensionality of B
  C = reshape(...
      repmat(reshape(A,prod(szA(1:end)),1),szB(1))*...
      reshape(B,szB(1),prod(szB(2:end))),...
      [szA(1:end), szB(1:end)]);
end

if dims == 1

  % Most important case (A last dimension matches B first)
  if szA(end) == szB(1)

    C = reshape(...
      reshape(A,prod(szA(1:end-1)),szA(end))*...
      reshape(B,szB(1),prod(szB(2:end))),...
      [szA(1:end-1),szB(2:end)]);

  % A last dimension doesn't match B first, but size(B,1) == 1
  elseif szB(1) == 1
    C = reshape(...
      reshape(A,prod(szA(1:end)),1)*...
      reshape(B,szB(1),prod(szB(2:end))),...
      [szA(1:end), szB(2:end)]);
    warning('The first dimension of B did not match the last dimension of A, so generalized outer product was assumed, extending A with singleton dimensions');
  else
    error('Gergo:mply:dimsDontMatch','Dimensions of the matrices do not match');
  end
end

% Take kind of tensor product with eliminating enough dimensions
if dims > 1
  %Make sure all the dimensions match
  for i1=1:dims
    if szA(end+1-i1) ~= szB(i1)
      error('Gergo:mply:dimsDontMatch','Dimensions of the matrices do not match');
    end
  end
  if dims == length(szA) % when getting rid of all A dimensions make sure there is a leading singleton
    A = reshape(A,[1,szA]);
    szA = size(A);
  end
  % Compute the lower dimensional product
  C = reshape(...
      reshape(A,prod(szA(1:(end-dims))),prod(szA((end-dims+1):end)))*...
      reshape(B,prod(szB(1:dims)),prod(szB((dims+1):end))),...
      [szA(1:(end-dims)),szB((dims+1):end)]);
end


% C = squeeze(C);
szC = size(C);
if sum(szC==1)>0
  szCsqueeze = szC(szC~=1);
  if numel(szCsqueeze)==1, szCsqueeze = [szCsqueeze,1]; end
  C = reshape(C,szCsqueeze);
end

% if length(szC)>2
%   if szC(1) == 1;
%     C = reshape(C,szC(2:end));
%   end
% end
