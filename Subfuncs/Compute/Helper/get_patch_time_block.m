function patch = get_patch_time_block( data, row,col, cutsize )
%GET_PATCH_TIME_BLOCK Summary of this function goes here
%   Detailed explanation goes here

szY = chomp_size(data.proc_stack,'Y');

d = floor(cutsize/2);


[ valid_inds, cuts ] = mat_boundary(szY(1:2), row-d:row+d, col-d:col+d);

patch = zeros(cutsize,cutsize,szY(3));

patch(1+cuts(1,1):end-cuts(1,2),1+cuts(2,1):end-cuts(2,2),:) = data.proc_stack.Y(valid_inds{1},valid_inds{2},:);


end

