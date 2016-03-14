function [W] = update_dict(datas,Hs,W,opts,k,type)
%UPDATE_DICT Summary of this function goes here
%   Detailed explanation goes here
% Use the mean image to find a reasonable backprojection of W to the
% original space, via the locations

if iscell(opts), opt=opts{1}; else opt = opts; end
  

switch opt.learn_decomp
  case 'LMSVD'
    patches = pick_patches(datas,Hs,opts,type);
    patches = flatten_patches(patches);
    [U, Sv] = lmsvd(patches,k);
  case 'MTF'
    %TODO Maybe sometime later
  case 'HOSVD'
    patches = pick_patches(datas,Hs,opts,type);
    patches = flatten_patches(patches);
    [U, Sv] = svd(patches,'econ');    
  otherwise %raise error
    error('CHOMP:learning:dict_update', 'Dictionary update option string (opt.learn_decomp) does not correspond to implemented options.');
    
end
    

if opt.fig >1
  figure(12);
  plot(diag(Sv)); title('singular values of HOSVD')
end

W(:,1:min(size(W,2),k)) = U(:,1:min(size(W,2),k));

% %Hacky translation step by Marius
% m = opt.m;
% d = floor(m/2);
% 
% W = reshape(W,m,m,[]);
% xs  = repmat(-d:d, m, 1);
% ys  = xs';
% 
% absW = abs(W);
% absW = absW/mean(absW(:));
% x0 =  mean2(mean(absW,3) .* xs);
% y0 =  mean2(mean(absW,3) .* ys);
% 
% xform = [1 0 0; 0 1 0; -x0 -y0 1];
% tform_translate = maketform('affine',xform);
% 
% for k = 1:size(W,3)
%     W(:,:,k) = imtransform(W(:,:,k), tform_translate,...
%         'XData', [1 m], 'YData',   [1 m]);
%     W(:,:,k) = W(:,:,k) - mean2(W(:,:,k));
%     if std2(W(:,:,k))~=0
%       W(:,:,k) = W(:,:,k)./std2(W(:,:,k));
%     end
% end
% 
% W = reshape(W,m^2,[]);

%Make sure all Ws are positive near the center
[~, mask] = transform_inds_circ(0,0,150,opt.m,max((opt.m-5)/2,1),0);
mask = logical(mask);
for k = 1:size(W,2)
  if sum(W(mask(:),k)) < sum(W(~mask(:),k))
    W(:,k) = -1.*W(:,k);
  end
end

function py=flatten_patches(patches)
  py = [];
  for i1 = 1:length(patches)
      out = [];
      patch = patches{i1};
      for mom = 1:opt.mom
        cur = patch{mom};
        cur = reshape(cur,opt.m^2,[]);
        cur = cur./size(cur,2); %normalize by dimensionality
        out = [out, cur];
      end
      py(:,end+1:end+size(out,2)) = out;
  end
end

end

