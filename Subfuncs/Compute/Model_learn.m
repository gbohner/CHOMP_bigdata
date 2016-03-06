function Model_learn( opt)
%MODEL_LEARN Summary of this function goes here
%   Detailed explanation goes here



load(get_path(opt), 'inp'); %could load opt as well, but it may have changed
inp.opt = struct_merge(inp.opt, opt);



if inp.opt.mask
  UserMask = inp.UserMask;
end

W  = Model_initialize(inp.opt);

tic;

%% Run the learning
for n = 1:inp.opt.niter 
    

    fprintf('Iteration %d/%d, inferring cell locations...\n', n, inp.opt.niter);
    %Compute convolution of Y with the filters as well as the "local Gram
    %matrices of filters to use in the matching pursuit step afterwards
    [WY, GW, WnormInv] = compute_filters(inp.data, W, inp.opt );
    
    %Infer the hidden variables (H - location, X - reconstruction weight, L - log likelihood gain)
    if ~inp.opt.mask
      [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, inp.opt);
    else
      [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, inp.opt, UserMask);
    end

    %Save results from current iteration
    results = chomp_results(inp.opt,W,H,X,L,inp.y,inp.y_orig,inp.V);
    save(get_path(inp.opt,'output_iter',n) ,'results')
    
    %Visualize results
    if inp.opt.fig >0
      update_visualize(results.y,results.H, ...
        reshape(results.W,results.opt.m,results.opt.m,size(results.W,ndims(results.W))),...
        results.opt,1);
    end

    if n < inp.opt.niter && inp.opt.learn
      %Update the dictionary (the W filters)
      fprintf('Iteration %d/%d, updating dictionary...\n', n, inp.opt.niter);
      [W] = update_dict(inp.data,results.H,results.W,inp.opt,n+2);
    end
    
   
   if rem(n,1)==0
        fprintf('Iteration %d/%d finished, elapsed time is %0.2f seconds\n', n, inp.opt.niter, toc)
    end
    

end

end

