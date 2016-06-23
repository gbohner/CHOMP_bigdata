function [ H, X, L] = extract_coefs( WY, GW, WnormInv, W, opt, varargin)
%EXTRACT_COV_COEFS Summary of this function goes here
%   Detailed explanation goes here

% opt.fig = 1; %TMP
if opt.fig > 1
  h_dl = figure(7);  
  h_dl2 = figure(8);  
  h_dl3 = figure(9);
  h_comb = figure('Units','Normalized','Position',[0.1,0.1,0.7,0.8]);
  if opt.fig >2
%     Video_dl = VideoWriter('likelihood_map.avi');
%     Video_dl.FrameRate = 2;
%     open(Video_dl);
%     Video_dl2 = VideoWriter('mean_score_map.avi');
%     Video_dl2.FrameRate = 2;
%     open(Video_dl2);
%     Video_dl3 = VideoWriter('var_score_map.avi');
%     Video_dl3.FrameRate = 2;
%     open(Video_dl3);
    Video_comb = VideoWriter('inference_video.avi');
    Video_comb.FrameRate = 2;
    open(Video_comb);
  end
end

Ntypes = opt.NSS;
m = opt.m;
sz = size(WY);

H = zeros(opt.cells_per_image,1); %Location (linearized index)
X = zeros(opt.cells_per_image, size(W,2) * opt.mom); % Basis function coefficients
L = zeros(opt.cells_per_image,1); % Likelihood gains


if opt.mask
  Mask = varargin{1};
else
  Mask = ones(size(WY,1),size(WY,2)); % Possible cell placements (no overlap / nearby cells);
end
Mask(1:opt.m,:) = 0; Mask(end-opt.m:end,:) = 0; Mask(:, 1:opt.m) = 0; Mask(:, end-opt.m:end) = 0;%Don't allow cells near edges
Mask = double(Mask);
%TODO: Maybe modify it such that mask is per object type


dL_mom = zeros([size(WY,1),size(WY,2),Ntypes,opt.mom]);
dL = zeros([size(WY,1),size(WY,2), Ntypes]); % delta log likelihood
xk = zeros([size(WY,1),size(WY,2), size(W,2), opt.mom]); % Coefficients for the mean image filter reconstruction


for j = 1:opt.cells_per_image
  
  % Update delta log-likelihoods
  if j == 1
    %Reset filter coefficients
    xk(:) = 0;

    %Update filter coefficients (MAP estimate)
    for map = 1:size(W,2)
      for mom = 1:opt.mom
        xk(:,:,:,mom) = xk(:,:,:,mom) + reshape(reshape(WY(:,:,map,mom),prod(sz(1:2)),1)*WnormInv(map,:,mom),[sz(1),sz(2),size(W,2)]);
      end
    end
  end
  
%     xk(xk<0) = 0; %TOTHINK
  
    %Compute delta log-likelihood blockwise
    for type = 1:Ntypes
      for mom = 1:opt.mom
        %Give relative weight to the moments based on how many elements
        %they involve
        dL_mom(:,:,type,mom) = - sum(WY(:,:,opt.Wblocks{type},mom) .* xk(:,:,opt.Wblocks{type},mom),3);
  %       if mom>=2
  %         dL_mom(:,:,mom) = dL_mom(:,:,mom)./abs(mean2(dL_mom(:,:,mom))); %normalize the moment-related discrepencies
  %       end
        dL_mom(:,:,type,mom) = reshape(zscore(reshape(dL_mom(:,:,type,mom),numel(dL_mom(:,:,type,mom)),1)),size(dL_mom(:,:,type,1)));
      end
%     dL = - sum(sum(WY .* xk,3),4); % Add contribution from each map and each moment
    
      dL(:,:,type) = sum(dL_mom(:,:,type,:),4); %linear sum of individual zscored differences %TODO GMM version of "zscoring jointly"
%     dL = dL_mom(:,:,4); % Make it kurtosis pursuit
    end

    % Find maximum decrease  
    [AbsMin, ind] = min( dL(:).*repmat(Mask(:),Ntypes,1) );
    [row, col, type] = ind2sub(size(dL),ind);
    
    %Store the values (before changing them below, omg...)
    H(j) = ind;
    for mom = 1:opt.mom
      X(j,(mom-1)*size(W,2)+opt.Wblocks{type}) = reshape(xk(row,col,(mom-1)*size(W,2)+opt.Wblocks{type}),1,[]); %TODO, check if leaving a bunch of zeros in X is reasonable, or rather have it just be size of opt.KSS and choose the respective block of W later
    end
    L(j,:) = AbsMin;
    
    %Check if there is not enough likelihood decrease anymore
    if AbsMin > 0
      break;
    end
  
    
    
    
  if opt.fig >1
%     set(0,'CurrentFigure',h_dl); imagesc(dL_mom(:,:,1)); colorbar; pause(0.05);
    set(0,'CurrentFigure',h_dl); imagesc(dL(:,:,1).*Mask); colorbar; axis square;  pause(0.05);
    set(0,'CurrentFigure',h_dl2); imagesc(dL_mom(:,:,1,min(1,size(dL_mom,4))).*Mask); colorbar; axis square; pause(0.05);
    set(0,'CurrentFigure',h_dl3); imagesc(dL_mom(:,:,1,min(2,size(dL_mom,4))).*Mask); colorbar; axis square; pause(0.05);
  
    if opt.fig >2
     load(get_path(opt,'output_iter',4),'model')
     set(0,'CurrentFigure',h_comb);
     subplot(2,2,1); imagesc(model.y_orig); colormap gray; axis square; axis off; title('Data');
     subplot(2,2,3); imagesc(dL(:,:,1).*Mask); colormap default; colorbar; axis square; axis off; title('Total cost change');  pause(0.05);
     subplot(2,2,2); imagesc(dL_mom(:,:,1,min(1,size(dL_mom,4))).*Mask); colorbar; axis square; axis off; title('r=1 cost change'); pause(0.05);
     subplot(2,2,4); imagesc(dL_mom(:,:,1,min(2,size(dL_mom,4))).*Mask); colorbar; axis square; axis off; title('r=2 cost change'); pause(0.05);
%    set(0,'CurrentFigure',h_dl3); imagesc(Mask(:,:,1)); colorbar; axis square; pause(0.05);
    end
  end
  
  
  
  
  %Affected local area
  % Size(,1) : number of rows, size(,2): number of columns
 [inds, cut] = mat_boundary(sz(1:2),row-m+1:row+m-1,col-m+1:col+m-1);
  
 

  % Compute the changes in WY and xk;
  for map = 1:size(W,2)
   for mom = 1:opt.mom
     for map2 = 1:size(W,2)
      WY(inds{1},inds{2},map2, mom) = WY(inds{1},inds{2},map2, mom) - ...
        GW{map,map2,mom}(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2))*xk(row,col,map+(mom-1)*size(W,2));
     end
   end
  end
 
  %Recompute the changed xk values
 
  xk(inds{1},inds{2},:,:) = 0; 
  for map = 1:size(W,2)
    for mom = 1:opt.mom
      xk(inds{1},inds{2},:,mom) = xk(inds{1},inds{2},:,mom) + reshape(reshape(WY(inds{1},inds{2},map,mom),numel(inds{1})*numel(inds{2}),1)*WnormInv(map,:,mom),[numel(inds{1}),numel(inds{2}),size(W,2)]);
    end
  end
 
%    figure(4); imagesc(WY(:,:,1)); colorbar; pause(0.05);  
 
  % Update the patch around the point found
%   Mask(max(row-3,1):min(row+3,end),max(col-3,1):min(col+3,end),type) = 0; % Make it impossible to put cells to close to eachother
%   Mask(max(row-1,1):min(row+1,end),max(col-1,1):min(col+1,end),:) = 0; % Make it impossible to put cells to close to eachother
    Mask(row,col,:) = 0; %Make it impossible to put a cell into the exact same location
  
if ~isempty(opt.spatial_push)
  [yinds, ycut] = mat_boundary(sz(1:2),row-opt.m:row+opt.m,col-opt.m:col+opt.m);  
  [gridx,gridy] = meshgrid((-opt.m+ycut(1,1)):(opt.m-ycut(1,2)),(-opt.m+ycut(2,1)):(opt.m-ycut(2,2))); %make sure to cut the corresponding dimensions using ycut
  gridx = gridx'; gridy = gridy'; %meshgrid(1:n,1:m) creates mxn matrices, need to transpose
  
  grid_dist = sqrt(gridx.^2+gridy.^2);
  grid_dist = opt.spatial_push(grid_dist); % Specified distance based function
  Mask(yinds{1},yinds{2},:) = Mask(yinds{1},yinds{2},:).*repmat(grid_dist,[1,1,size(Mask,3)]); % Make it impossible to put cells to close to eachother
end
  
 if opt.fig >2
%   writeVideo(Video_dl, getframe(h_dl));
%   writeVideo(Video_dl2, getframe(h_dl2));
%   writeVideo(Video_dl3, getframe(h_dl3));
  writeVideo(Video_comb, getframe(h_comb));
 end  
%   disp([num2str(j) ' cells found, current type: ' num2str(type)]);
end

  if opt.fig >2 
%     close(Video_dl);
%     close(Video_dl2);
%     close(Video_dl3);
    close(Video_comb);
  end
% close(Video_yres);
end

