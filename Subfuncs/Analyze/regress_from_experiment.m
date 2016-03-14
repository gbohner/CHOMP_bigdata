%%Loading stuff
% close all;
% %Load Dan's plotting tools if on neurofast
% addpath('/home/djoshea/code/rig1/analysis/djoshea')
% addpath('/home/djoshea/code/rig1/analysis/djoshea/utils')
% import Regress.plotTuningColorGuide

%Current time series
load(get_path(opt,'results'),'timeseries','ROIs','ROI_mask','opt');
timeseries_correct = timeseries';
%load('cur_time_series_rand', 'timeseries');
%timeseries_random = zscore(timeseries',[],1);
%Regression data
load('/neurotank/derived/Watkins/2016-02-12/2P_regression/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021/regressHandData.mat');

szY = chomp_size(inp.data.raw_stack,'Y');

%For timeseries just do the pixel-wise thing, see if I get the same results
%as Dan
%szY = size(inp.data.raw_stack,'Y');
%timeseries = reshape(inp.data.raw_stack.Y,prod(szY(1:2)),[])'; %pixelwise timeseries, each pixel is a column

% load('pixelwise_timeseries', 'timeseries');
% timeseries_correct = zscore(timeseries,[],1);

%% 
%Get rid of places where hand wasn't seen
timeseries_correct(((regressTable.handSeen~=1) + isnan(regressTable.handVelX) + isnan(regressTable.handVelY))>0,:) = [];
timeseries_random = timeseries_correct(randperm(size(timeseries_correct,1)),:);

timeseries_correct = zscore(timeseries_correct);
timeseries_random = zscore(timeseries_random);

%% Do the regression
%Getthe features from regressData we want to regress to
regressFrom = [regressTable.handVelX, regressTable.handVelY, ones(size(regressTable.handVelY))];
regressFrom(((regressTable.handSeen~=1) + isnan(regressTable.handVelX) + isnan(regressTable.handVelY))>0,:) = [];
regressFrom(:,1:2) = zscore(regressFrom(:,1:2));


reg_results = struct('b',[],'r',[],'r2',[]);

for cols =1:size(timeseries_correct,2)
  [reg_results(cols).b, ~, reg_results(cols).r, ~, stats]  = regress(timeseries_correct(:,cols), regressFrom);
  reg_results(cols).r2 = stats(1);
end


%% Do the regression to timeseries extracted from random ROIs 
reg_results_rand = struct('b',[],'r',[],'r2',[]);

for cols =1:size(timeseries_random,2)
  [reg_results_rand(cols).b, ~, reg_results_rand(cols).r, ~, stats] = regress(timeseries_random(:,cols), regressFrom);
  reg_results_rand(cols).r2 = stats(1);
end


%% Compare
disp('Residual sum squared normalized by std, first row: extracted time series, second row: random timeseries')
disp([[reg_results.r2]; [reg_results_rand.r2]])
figure; hist([[reg_results.r2]',[reg_results_rand.r2]'],32);

tmp = [reg_results.b];
tmp_rand = [reg_results_rand.b];
%Compute norm and direction of coefficients for each cell
%reg_results.

vel_pred = sqrt(sum(tmp.^2));
vel_pred_rand = sqrt(sum(tmp_rand.^2));
dir_pred = atan2(tmp(2,:),tmp(1,:));
dir_pred_rand = atan2(tmp_rand(2,:),tmp_rand(1,:));

disp('Prediction of activity by ROI, hand speed for original and randomized timeseries')
disp([vel_pred; vel_pred_rand])
mean([vel_pred; vel_pred_rand],2)

%% Color the ROIs according to regression parameters
load('cur_time_series', 'ROIs');
%pr = load('pixelwise_reg');
%tmp = [pr.reg_results.b];
%vel_pred_pix = sqrt(sum(tmp.^2));


out_im = zeros([szY(1:2),3]);
out_im_rand = zeros([szY(1:2),3]);
for i1 = 1:numel(ROIs)
  hsvList = [(dir_pred(i1)'+pi)/(2*pi), repmat(0.5, 1, 1), mat2gray(vel_pred(i1), [0,1])]; %[min(vel_pred_pix(:)), max(vel_pred_pix(:))])];
  rgbList = hsv2rgb(hsvList);
  hsvList_rand = [(dir_pred_rand(i1)'+pi)/(2*pi), repmat(0.5, 1, 1), mat2gray(vel_pred_rand(i1), [0,1])]; %[min(vel_pred_pix(:)), max(vel_pred_pix(:))])];
  rgbList_rand = hsv2rgb(hsvList_rand);
  col = ROIs{i1}.col;
  row = ROIs{i1}.row;
  colored_ROI = reshape(reshape(ROIs{i1}.mask,1,[])'*rgbList,[size(ROIs{i1}.mask),3]);
  colored_ROI_rand = reshape(reshape(ROIs{i1}.mask,1,[])'*rgbList_rand,[size(ROIs{i1}.mask),3]);
 
  cutsize = size(ROIs{i1}.mask,1);
  [inds, cut] = mat_boundary(szY(1:2),row-floor(cutsize/2):row+floor(cutsize/2),col-floor(cutsize/2):col+floor(cutsize/2));
  out_im(inds{1},inds{2},:) = colored_ROI(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2),:);
  out_im_rand(inds{1},inds{2},:) =  colored_ROI_rand(1+cut(1,1):end-cut(1,2),1+cut(2,1):end-cut(2,2),:);
end

% figure;
% subplottight(2, 1, 1);
% imagesc(imread('pixelwise_image.png')); 
% subplottight(2, 1, 2);
% imagesc(mat2gray(y_orig)); colormap gray
% axes('Position', [0.85 0.85 0.15 0.15]);
% plotTuningColorGuide();
% print(gcf,'pixelwise_map.png','-dpng')

figure;
subplottight(2, 1, 2);
imagesc(out_im_rand);
subplottight(2, 1, 1);
imagesc(out_im);
axes('Position', [0.85 0.85 0.15 0.15]);
plotTuningColorGuide();
%print(gcf,'ROI_colors_map.png','-dpng')




%% Color image according to pixelwise regression results
load('pixelwise_reg');
tmp = [reg_results.b];
vel_pred_pix = sqrt(sum(tmp.^2));
dir_pred_pix = atan2(tmp(2,:),tmp(1,:));

hsvList = [(dir_pred_pix'+pi)/(2*pi), repmat(0.5, numel(dir_pred_pix), 1), mat2gray(vel_pred_pix')];
rgbList = hsv2rgb(hsvList);
rgb = reshape(rgbList, szY(1), szY(2), 3);
figure; imagesc(rgb); axis image
%imwrite(rgb, 'pixelwise_image.png');
