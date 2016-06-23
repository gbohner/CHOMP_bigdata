function [out, num_cells, col_count] = pick_patches( datas, Hs, opts,type, varargin)
%PICK_PATCHES Returns patches to learn from, Y is the data tensor, H are the spatial locations, type is cell type number

%Make sure all the input are cell arrays (of possibly multiple datasets)
if ~iscell(datas), datas = {datas}; end
if ~iscell(Hs), Hs = {Hs}; end
if ~iscell(opts), opts = {opts}; end

do_cov = 2; 
% 0 - return all n order patches, 
% 1 - return covariance matrix build from n-order patches, 
% 2- return covariance matrix built from raw patches
if nargin > 4 %Check if we want to just output the covariance matrix instantly
  do_cov = varargin{1};
end

py = cell(numel(datas),1);
num_cells = zeros(numel(datas,1));

parfor c1 = 1:numel(datas)
  data = datas{c1};
  H = Hs{c1};
  opt = opts{c1};
  
  szY = chomp_size(data.proc_stack,'Y');
  
  %Remove entries from H that are the wrong type
  h1 = 1;
  while h1<=numel(H)
    [~,~,cur_type] = ind2sub([szY(1:2) opt.NSS],H(h1));
    if cur_type ~= type, H(h1) = []; else h1 = h1+1; end
  end
  
  
  if do_cov
    py{c1} = struct('mat',zeros(opt.m^2), 'count', 0);
  else
    py{c1} = cell(numel(H),1);
  end
  
  if ~isempty(H)
    patches = get_patch(data.proc_stack, opt, H);
  else
    continue;
  end
  
  num_cells(c1,1) = numel(H);
  
  if do_cov == 2
    %Just build the covariance matrix out of the raw patch samples
    patches = reshape(patches,opt.m^2,[]);
    patches = bsxfun(@minus, patches, mean(patches,2));
    py{c1}.count = size(patches,2);
    py{c1}.mat = (patches * patches');
  else
    %Process the individual patches to get n-order estimates
    for h1 = 1:numel(H)
      %disp(h1);
      curpy = get_n_order_patch(patches(:,:,:,h1), opt, szY);
      %Also weigth every moment tensor according to the number of independent
      %elements (patchsize multichoose mom) over total number of elements
      for mom1 = 1:opt.mom
        curpy{mom1} = curpy{mom1} .* (nchoosek(opt.m^2+mom1-1,mom1)./((opt.m^2).^mom1));
      end
      if do_cov
        %Just store the resulting covariance matrix
        for mom1 = 1:opt.mom
          py{c1}.mat = py{c1}.mat + ...
            reshape(curpy{mom1},size(curpy{mom1},1),[])*reshape(curpy{mom1},size(curpy{mom1},1),[])'; %reshape to flat for "HOSVD"
          py{c1}.count = py{c1}.count + size(reshape(curpy{mom1},size(curpy{mom1},1),[]),2); %number of columns
        end
      else
        %Store all the individual, weighted and flattened vectors
        py{c1}{h1} = curpy;
      end
    end
  end


end

opt = opts{1};

if do_cov
  out = zeros(opt.m^2, opt.m^2);
  col_count = 0;
  for c1 = 1:numel(py);
    %Combine all the info from all datasets, with numerical stability
    weigth = (py{c1}.count - 1);
    out = out + py{c1}.mat./weigth; %number of samples
    col_count = col_count + py{c1}.count;
  end
else
  %Concatanate the results into a
  out={};
  for c1 = 1:numel(py)
    out(end+1:end+numel(py{c1})) = py{c1};
  end

  out = flatten_patches(out, opt); % opt.m^2 x (location*(opt.m^2)^opt.mom)  - very flat matrix
  col_count = size(out,2);
end

num_cells = sum(num_cells);


function py=flatten_patches(patches,opt)
  py = [];
  for i1 = 1:length(patches)
      out1 = [];
      patch = patches{i1};
      for mom = 1:opt.mom
        cur = patch{mom};
        cur = reshape(cur,opt.m^2,[]);
        cur = cur./size(cur,2); %normalize by dimensionality
        out1 = [out1, cur];
      end
      py(:,end+1:end+size(out1,2)) = out1;
  end
end

end