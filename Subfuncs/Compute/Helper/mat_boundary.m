function [ valid_ind, cuts ] = mat_boundary(sz, varargin)
%BOUNDARY takes a matrix size and a range of indices as input and outputs a
%range of valid indices, as well as how much cut was done in each dimension
%   2D for now

cuts = zeros(length(sz),2);
valid_ind = cell(length(sz),1);

for i1 = 1:length(sz)
  [valid_ind{i1}, cuts(i1,1), cuts(i1,2)] = b_1d(sz(i1),varargin{i1});
end
  
  
  function [valid_ind1, cmin, cmax] =  b_1d(n,ind1)
    cmin = sum(ind1<1); cmax = sum(ind1>n); %How much it is out of boundary
    valid_ind1 = ind1(ind1>=1 & ind1 <= n);
  end


end

