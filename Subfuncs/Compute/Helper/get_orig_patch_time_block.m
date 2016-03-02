function patch = get_orig_patch_time_block( data, row,col, opt )
%GET_PATCH_TIME_BLOCK Summary of this function goes here
%   Detailed explanation goes here

szRaw = chomp_size(data.raw_stack,'Y');

%Rescale params back to original size
row = round(row./opt.spatial_scale);
col = round(col./opt.spatial_scale);
cutsize = round(opt.m./opt.spatial_scale);
cutsize = cutsize + (1-mod(cutsize,2)); %make sure it is an odd number

d = floor(cutsize/2);


[ valid_inds, cuts ] = mat_boundary(szRaw(1:2), row-d:row+d, col-d:col+d);

patch = zeros(cutsize,cutsize,szRaw(3));

patch(1+cuts(1,1):end-cuts(1,2),1+cuts(2,1):end-cuts(2,2),:) = ...
  data.raw_stack.Y(valid_inds{1},valid_inds{2},:);

%Using the virtual stack of the raw data

%{
% Using file reading all the time
data.raw.ReadFcn =  @(filename)imread(filename,'PixelRegion', ...
  {[valid_inds{1}(1), valid_inds{1}(end)], [valid_inds{2}(1), valid_inds{2}(end)]});

tmp_patch = readall(data.raw);
patch(1+cuts(1,1):end-cuts(1,2),1+cuts(2,1):end-cuts(2,2),:) = shiftdim(reshape(cell2mat(tmp_patch),[szRaw(3),numel(valid_inds{1}), numel(valid_inds{2})]),1);
%patch(1+cuts(1,1):end-cuts(1,2),1+cuts(2,1):end-cuts(2,2),:) = shiftdim(reshape(cell2mat(readall(data.raw)),[szRaw(3),numel(valid_inds{1}), numel(valid_inds{2})]),1);

data.raw.ReadFcn = @(filename)imread(filename);
%}
end

