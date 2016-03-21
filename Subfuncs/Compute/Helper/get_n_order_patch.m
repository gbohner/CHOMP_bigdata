function patch_out = get_n_order_patch(patch_block, opt, szY)    
  patch_block = reshape(patch_block,opt.m^2,[]); %patchsize * T
  patch_out = cell(1,opt.mom);
  tmp = cell(1,opt.mom);
  for t1 = 1:size(patch_block,2)
    for mom = 1:opt.mom
      if mom==1, tmp{mom} = patch_block(:,t1); else
        tmp{mom} = mply(tmp{mom-1},patch_block(:,t1)',0); %take outer product to compute next raw moment
      end
      if t1 == 1, patch_out{mom} = zeros(size(tmp{mom})); end; %initialize as 0s
      patch_out{mom} = patch_out{mom} + tmp{mom}./szY(end); %Add the patch's momth moment tensor divided by total time points
    end
  end
  
  %Convert the raw moment tensors into cumulant tensors
  patch_out = raw2cum_multivariate( patch_out );
end

