%% Loading
%Load the input file we wanna use (mostly manual for now)
%load(get_path(opt),'inp'); 
setenv('CHOMP_ROOT_FOLDER','/mnt/stanford')

addpath('/mnt/stanford/home/djoshea/code/rig1/analysis/djoshea')
addpath('/mnt/stanford/home/djoshea/code/rig1/analysis/djoshea/utils')
import Regress.plotTuningColorGuide

load('/mnt/stanford/neurotank/derived/gbohner/input/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021_Cycle00001_Ch2_000001.ome_20160310T134413.mat')
data = inp.data;
load(get_path(inp.opt,'output_iter',inp.opt.niter),'model');
[H, W, X, y_orig, y] = model.get_fields( 'H', 'W', 'X', 'y_orig','y');

opt=inp.opt;



%% Getting ROIs
close all;
getRandom = 0;
update_visualize( y,H,reshape(W,opt.m,opt.m,size(W,2)),opt,1,1);
%opt.ROI_type = 'mean_origsize';
%opt.ROI_type = 'quantile_origsize';
%opt.ROI_type = 'quantile';
opt.ROI_type = 'quantile_dynamic_origsize';
opt.ROI_params = 0.6;
%opt.ROI_params = 0.7;
if getRandom
  [ROI_mask, ROIs] = getROIs(opt, min(30,numel(H)),1); opt.fig = 2;
else
  [ROI_mask, ROIs] = getROIs(opt, min(30,numel(H)),0); opt.fig = 2;
end
figure(5);
% subplottight(2, 1, 2);
h_rois = imagesc(y_orig); colormap gray; axis image;
B = bwboundaries(ROI_mask);
hold on;
visboundaries(B)
for i1 = 1:numel(ROIs)
  text(ROIs{i1}.col, ROIs{i1}.row, num2str(i1), 'Color','r','FontSize',20,'FontWeight','bold');
end
  set(gca, 'XTick', []);
   set(gca, 'YTick', []);
%    
% subplottight(2, 1, 1);
% imagesc(out_im); axis image
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% axes('Position', [0.85 0.85 0.15 0.15]);
% plotTuningColorGuide();
 
pause(0.3);

if getRandom
  print(gcf,'cur_ROIs_rand.eps','-depsc2')
else
  print(gcf,'cur_ROIs.eps','-depsc2')
  print(gcf,'cur_ROIs.png','-dpng')
end


%% Getting timeseries from ROIs
szY = chomp_size(data.raw_stack, 'Y');
%timeseries = zeros(numel(H), szY(3));

patches = get_patch(data.raw_stack,opt,H,1:szY(3));

for i1 = 1:numel(H)
  timeseries(i1,:) = mply(ROIs{i1}.mask, patches(:,:,:,i1),2);
end

if getRandom
  save('cur_time_series_rand','timeseries','ROIs','ROI_mask');
else
  save('cur_time_series','timeseries','ROIs','ROI_mask');
end

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