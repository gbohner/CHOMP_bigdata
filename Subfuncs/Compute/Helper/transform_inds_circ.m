function [newinds, mask] = transform_inds_circ( c1,c2,size1, fs, mask1r, mask2r )
%TRANSFORM_INDS Summary of this function goes here
%   c1-c2 top left corner coordinate, size1 is size(matrix,1), fs is
%   filtersize, the number of pixels

  newinds = zeros(fs.^2,1);
  newinds(1) = transform_single_ind(c1,c2,size1);
  newinds(1:fs) = (newinds(1):(newinds(1)+fs-1))';
  newinds = repmat(newinds(1:fs),fs,1) + reshape(repmat((0:(fs-1))*size1,fs,1), fs.^2,1);
  
  mask = padarray(fspecial('disk',mask1r),[floor((fs-(mask1r*2+1))/2), floor((fs-(mask1r*2+1))/2)]);
  if mask2r>0
    mask2 = padarray(fspecial('disk',mask2r),[floor((fs-(mask2r*2+1))/2), floor((fs-(mask2r*2+1))/2)]);
    mask(mask2>0) = 0;
  end
  mask(mask>0) = 1;
  
  
  function out = transform_single_ind(c1,c2,size1)
    out = (c2-1) * size1 + c1;
  end

end

