%% Loading
load(get_path(opt),'inp');
data = inp.data;
load(get_path(inp.opt,'output_iter',inp.opt.niter),'results');
[H, W, X, y_orig, y] = results.get_fields( 'H', 'W', 'X', 'y_orig','y');

opt=inp.opt;

%% Getting ROIs
close all;
update_visualize( y,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1,1);
%opt.ROI_type = 'mean_origsize';
%opt.ROI_type = 'quantile_origsize';
%opt.ROI_type = 'quantile';
opt.ROI_type = 'quantile_dynamic_origsize';
opt.ROI_params = 0.6;
%opt.ROI_params = 0.7;
[ROI_mask, ROIs] = getROIs(opt, min(30,numel(H))); opt.fig = 2;
figure(5); imagesc(y_orig); colormap gray; axis image;
B = bwboundaries(ROI_mask);
hold on;
visboundaries(B)
pause(0.3);

%% Getting timeseries from ROIs
szY = chomp_size(data.proc_stack, 'Y');
timeseries = zeros(numel(H), szY(3));

for i1 = 1:numel(H)
    [row, col] = ind2sub(size(y),H(i1));
    %TODO something is being weird with reading patches from the original
    %image stack, prob fixed by matlab virtual stacks of raw data.
    patch = get_orig_patch_time_block( data, row,col, opt ); %this is the correct version (cause of imread being weird)
    %patch = get_patch_time_block( data, row,col, opt.m ); %this is the correct version (cause of imread being weird)
    timeseries(i1,:) = ROIs{i1}.mask(:)'*reshape(patch,size(patch,1)*size(patch,2),[]);
%     figure(6); imagesc(mean(patch,3))
%     figure(7); imagesc(ROIs{i1}.mask)
%     pause
    %pause(0.5);
end

save('cur_time_series','timeseries','ROIs','ROI_mask');

%% Just plotting

figure; 
to_plot = [1:20];%[15:20];%[10:15]+30;
v = max(std(timeseries(to_plot,:),1))*2;
for i1 = to_plot
  plot(timeseries(i1,:) + numel(to_plot)*v - i1*v, 'LineWidth', 2); hold on;
  set(gca,'YTick',[])
  xlabel('Frame')
end


%% Plotting the PCs of the timeseries
[coeff,score,latent] = pca(timeseries,'Algorithm', 'svd', 'Centered',true);
%[coeff,score,latent] = pca(timeseries);
[u,s,v]=  svd(bsxfun(@times, bsxfun(@minus, timeseries, mean(timeseries,1)),1./std(timeseries,[],1))');
%[u,s,v]=  svd(bsxfun(@minus, timeseries, mean(timeseries,1))');
[u,s,v] = svd(timeseries');
%[u,s,v1] = eig( zscore(timeseries)' * zscore(timeseries));
coeff = coeff';
figure;
to_plot = 1:5;
v = max(std(coeff(to_plot,:),1))*3;
for i1 = to_plot
  plot(coeff(i1,:) + numel(to_plot)*v - i1*v); hold on;
  set(gca,'YTick',[])
  xlabel('Frame')
end