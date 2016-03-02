function Model_learn( opt)
%MODEL_LEARN Summary of this function goes here
%   Detailed explanation goes here



load(get_path(opt), 'y', 'y_orig', 'data'); %could load opt as well, but it may have changed

if opt.mask
  load(get_path(opt), 'UserMask');
end

W  = Model_initialize(opt);

tic;

%% Run the learning
for n = 1:opt.niter 
    

    fprintf('Iteration %d/%d, inferring cell locations...\n', n, opt.niter);
    %Compute convolution of Y with the filters as well as the "local Gram
    %matrices of filters to use in the matching pursuit step afterwards
    [WY, GW, WnormInv] = compute_filters(data, W, opt );
    
    %Infer the hidden variables (H - location, X - reconstruction weight, L - log likelihood gain)
    if ~opt.mask
      [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, opt);
    else
      [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, opt, UserMask);
    end

    %Save results from current iteration
    save(get_path(opt,'output_iter',n) ,'y','y_orig','H','X','L','W','opt')% 
    
    %Visualize results
    if opt.fig >0
      update_visualize( y,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1);
    end

    if n < opt.niter && opt.learn
      %Update the dictionary (the W filters)
      fprintf('Iteration %d/%d, updating dictionary...\n', n, opt.niter);
      [W] = update_dict(data,H,W,opt,n+2);
    end
    
   
   if rem(n,1)==0
        fprintf('Iteration %d/%d finished, elapsed time is %0.2f seconds\n', n, opt.niter, toc)
    end
    

end

end

