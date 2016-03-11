function patch_out = get_n_order_patch(patch_block, opt, szY)    
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

