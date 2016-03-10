function [py] = pick_patches( data, H, opt, type)
%PICK_PATCHES Returns patches to learn from, Y is the data tensor, H are the spatial locations, type is cell type number

szY = chomp_size(data.proc_stack,'Y');

i1 = 1;
while i1<=size(H,1)
  [row,col,t] = ind2sub(szY(1:2),H(i1));
  if (t~=type)
    H(i1) = [];
  else
    i1 = i1+1;
  end
end

py = {};

patches = get_patch(data.proc_stack, opt, H);

%Pick the patches
for h1 = 1:size(H,1)
  curpy = get_n_order_patch(patches(:,:,:,h1));
  py{end+1} = curpy;
end

  function patch_out = get_n_order_patch(patch_block)    
    patch = reshape(num2cell(reshape(patch_block,opt.m^2,[]),1),szY(end),1);
    patch = repmat(patch,[1,opt.mom]); %T * moments cell
    patch_out = cell(1,opt.mom);
    for mom = 1:opt.mom
      for t1 = 1:length(patch);
        if mom>1, patch{t1,mom} = mply(patch{t1,mom-1},patch{t1,1}',0); end
        if t1 == 1, patch_out{mom} = zeros(size(patch{t1,mom})); end; %initialize as 0s
        patch_out{mom} = patch_out{mom} + patch{t1,mom}./szY(end); %Add the patch's momth moment tensor divided by total time points
      end
    end
  end
            

  

end

