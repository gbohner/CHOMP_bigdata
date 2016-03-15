function [WY, GW, WnormInv] = compute_filters(data, W, opt )
%COMPUTE_FILTERS Computes the correlation between filter and data, plus the
% original MAP coefficients
%   Detailed explanation goes here

%Load the GPT (we dont really want to keep it memory all the time
if ~exist(get_path(opt, 'precomputed'),'file')
  GPT = precompute_shift_tensors(opt);
else
  load(get_path(opt, 'precomputed'), 'GPT');
end

szY = chomp_size(data.proc_stack,'Y'); %size of the data tensor

WY = zeros([szY(1:2),opt.NSS*opt.KS,opt.mom]); %big matrix storing the vector resulting of filtering for each location, filter and moment

% Compute the convolutions with the filters
for filt = 1:size(W,2); %Different filters
  Wcur = W(:,filt);
  Wcurc = Wcur(:);
  Wcurc = Wcurc./norm(Wcurc+1e-6); % make sure it has norm of 1.
  Wconv = reshape(Wcurc,opt.m,opt.m);
  
  for t1 = 1:szY(end) %Can be done parallelly or on GPU 
    conv_result = conv2(data.proc_stack.Y(:,:,t1),Wconv,'same');
    for mom = 1:opt.mom  %Get raw moments of the projected time course at each possible cell location %TODO - this might be wrong, because it assumes equal weighting??? but it's filter by filter, so maybe the linear combination of filters is still linear in the higher moments %TOTHINK Nah it seems correct
      WY(:,:,filt,mom) = WY(:,:,filt,mom) + conv_result.^mom./szY(end);
    end
  end
end

%Convert the diagonal raw  moment estimates into diagonal cumulant
%estimates %TODO: first compute the full non-diag raw moments, convert to
%non-diag cumulants to estimate error coming from using diag versions
for type = 1:opt.NSS
  WY(:,:,opt.Wblocks{type},:) = raw2cum(WY(:,:,opt.Wblocks{type},:),4);
end


WnormInv = zeros(size(W,2),size(W,2),opt.mom); % Inverse Interaction between basis functions

%WnormInv should be in feature space (when doing random projection version)
for type = 1:opt.NSS
  for filt1 = 1:opt.KS
    for filt2 = 1:opt.KS
      Wcur1 = W(:,(type-1)*opt.KS+filt1);
      Wcur2 = W(:,(type-1)*opt.KS+filt2);
      for mom = 1:opt.mom
        if mom>1, Wcur1 = mply(Wcur1, W(:,(type-1)*opt.KS+filt1)',0); end
        if mom>1, Wcur2 = mply(Wcur2, W(:,(type-1)*opt.KS+filt2)',0); end
        Wcur1 = Wcur1./norm(Wcur1(:)+1e-6); % make sure it has norm of 1.
        Wcur2 = Wcur2./norm(Wcur2(:)+1e-6); % make sure it has norm of 1.

        WnormInv((type-1)*opt.KS+filt1,(type-1)*opt.KS+filt2,mom) = Wcur1(:)'*Wcur2(:); 
      end
    end
  end
end

%TODO: Something is wrong with reconstructions

%Invert Wnorm blockwise / obj type
for mom = 1:opt.mom
  for type = 1:opt.NSS
    WnormInv(opt.Wblocks{type},opt.Wblocks{type},mom) = ...
      inv(WnormInv(opt.Wblocks{type},opt.Wblocks{type},mom)+eye(numel(opt.Wblocks{type}))); % Regularized
  end 
end


% Use Worig to compute the matching pursuit step in the original space
GW = cell(size(W,2), size(W,2), opt.mom); %Each filter combination at each moment


%Each cell is going to be a cell of (2*m-1)^2 shifts and at each shift and
%each moment we'll have a vector of features^moment to describe how much
%the corresponding WY entry is modified if we set the coeffecient of active
%filt1 at moment mom to 1.

Worig = W;

for filt1 = 1:size(W,2)
  for filt2 = 1:size(W,2)
    Wcur1 = Worig(:,filt1);
    Wcur2 = Worig(:,filt2);
    for mom = 1:opt.mom
      if mom>1, Wcur1 = mply(Wcur1, Worig(:,filt1)',0); end
      if mom>1, Wcur2 = mply(Wcur2, Worig(:,filt2)',0); end

      %TODO flip all dimensions of the second filter, such that convolution
      %gives you nd correlation instead
%       Wcur2r = Wcur2;
%       Wcur2r = flipdim_all(Wcur2r);
      Wcur1 = Wcur1./norm(Wcur1(:)+1e-6); % make sure it has norm of 1.
      Wcur2 = Wcur2./norm(Wcur2(:)+1e-6); % make sure it has norm of 1.
      Wcur1c = Wcur1(:);
      Wcur2c = Wcur2(:);

      GW{filt1,filt2,mom} = zeros(2*opt.m-1, 2*opt.m-1);
      for s1 = 1:(2*opt.m-1)
        for s2 = 1:(2*opt.m-1)
          GW{filt1,filt2,mom}(s1,s2) = Wcur2c(GPT{s1,s2,mom}(:,2))'* Wcur1c(GPT{s1,s2,mom}(:,1)); %compute the shifted effect in original space via the shift tensors GPT. Because the Worigs were computed to correspond to the best inverse of the Ws
        end
      end
    end
  end
end

  
clearvars -except WY GW WnormInv

end

