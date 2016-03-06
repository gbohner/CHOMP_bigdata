close all;
clear all;

cd(fileparts(mfilename('fullpath')));

cur_opt = struct(...
    'data_path', '~/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif', ...
    'data_type', 'frames_virtual', ...
    'init_model','filled', ...
    'niter', 4, ...
    'm', 17, ...
    'mom', 2, ...
    'fig',1, ...
    'spatial_scale',1,...
    'time_scale',1,...
    'mask', 0, ...
    'smooth_filter_mean',2, ...
    'smooth_filter_var',2, ...
    'cells_per_image', 50, ...
    'KS', 4 ...
  );

[ROI_mask, ROIs] = chomp(cur_opt);

% load('../data_for_gergo/S1-T52844_masks.mat', 'donut')

% figure(2); imshowpair(ROI_mask, donut>0);