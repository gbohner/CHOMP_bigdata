function [W,  Worig]  = Model_initialize( opt )
%MODEL_INITIALIZE Initializes the basis function matrices
%   Detailed explanation goes here
  
  W = zeros(opt.m^2, opt.NSS*opt.KS); % initialize basis functions
    % Initialize circles with reasonable size
  
  for type = 1:opt.NSS
    init_model = opt.init_model{type};
    switch init_model
      case 'filled'
        [~, mask] = transform_inds_circ(0,0,150,opt.m,(opt.m-1)/2,0); % . , . , ., filter size, circle outer radius, inner hole radius
        W(:,opt.Wblocks{type}(1)) = mask(:);
      case 'supervised'
        %learn the best basis functions from varargin{2}, which should be a
        %collection of examples or alternative a set of locations within
        %the dataset (in which case just call an update_dict with the
        %preset H)
      case 'pointlike'
        % Initialize the to a dot/small circle
        [~, mask] = transform_inds_circ(0,0,150,opt.m,min((opt.m-1)/2,3),0);
        W(:,opt.Wblocks{type}(1)) = mask(:);
      case 'donut'
        [~, mask] = transform_inds_circ(0,0,150,opt.m,(opt.m-1)/2,(opt.m-5)/2); % . , . , ., filter size, circle outer radius, inner hole radius
        W(:,opt.Wblocks{type}(1)) = mask(:);
      case 'given'
        W(:,opt.Wblocks{type}) = opt.init_W(:,opt.Wblocks{type});
      otherwise
        error('CHOMP:learning:initialize', 'Model initialization option string (opt.init_model) does not correspond to implemented options.');
    end
  end
    
%      [~, mask] = transform_inds_circ(0,0,150,opt.m,(opt.m-3)/2,0); % . , . , ., filter size, circle outer radius, inner hole radius
%     W(:,2) = mask(:);
%      [~, mask] = transform_inds_circ(0,0,150,opt.m,(opt.m-5)/2,0); % . , . , ., filter size, circle outer radius, inner hole radius
%     W(:,3) = mask(:);
  
  
  Worig = W;
  
  %Project the assumed initial basis function into feature space
%   if opt.rand_proj
%     W = opt.P*W;
%   end
  
  % Make sure that basis function column norms are ~1. 
  W = bsxfun(@times, W, 1./sum(W.^2+1e-6,1));
  
 
  
end

