function out = get_patch( stack, opt, H, varargin )
%GET_PATCH Gets patch(es) from the data at given location at time
%   get_patch(data,opt,H) - returns an array of patch time blocks
%   around locations H, size specified as in opt.m, out = row x col x t x
%   cell_num
%   get_patch(data,opt,H, t) - returns a cell array of patches blocks
%   around locations H at time t, size specified as in opt.m

szY = chomp_size(stack,'Y');

if nargin==3
  times = 1:szY(end);
else
  times = varargin{1};
end

szPatch = [opt.m, opt.m];

%Initialize output array
out = cell(numel(times),1);
[out{:}]= deal(zeros([numel(H),szPatch]));


for t = times %over frames required
  frame_t = stack.Y(:,:,t);
  for i1 = 1:numel(H)
    [row, col] = ind2sub(szY,H(i1));
    [ valid_inds, cuts ] = mat_boundary(szY(1:2), row-floor(szPatch(1)/2):row+floor(szPatch(1)/2), col-floor(szPatch(1)/2):col+floor(szPatch(1)/2));
    out{t}(i1,1+cuts(1,1):end-cuts(1,2),1+cuts(2,1):end-cuts(2,2)) = frame_t(valid_inds{1},valid_inds{2});
  end
end

out = reshape(cell2mat(out),[numel(H),numel(times),szPatch]);
out = permute(out, [3,4,2,1]); % Change to row x col x t x cell_num
end

