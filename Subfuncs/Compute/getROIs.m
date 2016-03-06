function [ ROI_image, ROIs ] = getROIs( opt, varargin )
%GETROIS Summary of this function goes here
%   Detailed explanation goes here

load(get_path(opt, 'output_iter', opt.niter), 'results');
 [H, W, X, y_orig, y] = results.get_fields( 'H', 'W', 'X', 'y_orig','y');

%update_visualize( y_orig,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1);

if nargin>1
  num_reconst = varargin{1};
else
  num_reconst = length(H);
end

sz = size(y);

ROIs = cell(num_reconst,1);
switch opt.ROI_type
    case 'quantile'
      ROI_image = zeros(sz(1:2));
    case 'quantile_dynamic'
      ROI_image = zeros(sz(1:2));
    case 'quantile_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
    case 'quantile_dynamic_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
    case 'mean_origsize'
      szRaw = size(y_orig);
      ROI_image = zeros(szRaw(1:2));
end

for i1 = 1:num_reconst
  [row,col,t] = ind2sub(sz,H(i1));
  reconst = reshape(W(:,1:opt.KS)*(X(i1, 1:opt.KS)'), opt.m, opt.m);
  reconst = imrotate(reconst, 180);
  
  switch opt.ROI_type
    case 'quantile'
      reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      reconst = bwconvhull(reconst);
      [inds, cut] = mat_boundary(sz(1:2),row-floor(opt.m/2):row+floor(opt.m/2),col-floor(opt.m/2):col+floor(opt.m/2));
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'quantile_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      reconst = imdilate(reconst,strel('rectangle',[3,3]));
      reconst = imdilate(reconst,strel('rectangle',[3,3]));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      reconst = imerode(reconst,strel('rectangle',[3,3]));
      %reconst = bwconvhull(reconst); %looks better for visualizations
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'mean_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      %reconst = reconst > quantile(reconst(:),opt.ROI_params(1));
      %reconst = bwconvhull(reconst); %looks better for visualizations
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
%       %get the covariance between pixels of reconstruction and original
%       %mean
%       reconst = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2)) .* y_orig(inds{1},inds{2});
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    case 'quantile_dynamic_origsize'
      cutsize = round(opt.m./opt.spatial_scale);
      cutsize = cutsize + (1-mod(cutsize,2));
      reconst = imresize(reconst,[cutsize,cutsize]);
      row = round(row./opt.spatial_scale);
      col = round(col./opt.spatial_scale);
      [inds, cut] = mat_boundary(szRaw(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
      tmp = opt.ROI_params(1)*(max(reshape(y_orig(inds{1},inds{2}),1,[])')-min(reshape(y_orig(inds{1},inds{2}),1,[])'))+min(reshape(y_orig(inds{1},inds{2}),1,[])');
      tmp = sum(sum(y_orig(inds{1},inds{2})>tmp))./numel(y_orig(inds{1},inds{2})); %ratio of pixels to pick
      [~, tmp_ind] = sort(reconst(:),'descend');
      reconst(:) = 0;
      reconst(tmp_ind(1:floor(tmp*numel(reconst))))=1;
      ROI_image(inds{1},inds{2}) = reconst(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2));
    otherwise
      error('CHOMP:roi:method',  'Region of interest option string (opt.ROI_type) does not correspond to implemented options.')
  end
  
  %Store results in cells now, later add option for json output %TODO
  ROIs{i1} = struct('col', col, 'row', row, 'mask', reconst);
  
  
  
  
end

  if opt.fig >1
    figure(2);
    imshow(ROI_image); axis square;
  end

end

