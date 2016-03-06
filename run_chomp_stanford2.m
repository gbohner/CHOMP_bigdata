close all;
clear all;

cd(fileparts(mfilename('fullpath')));
addpath('./Classes')

opt = chomp_options(...
    'root_folder', '/mnt/stanford', ...
    'input_folder', '/neurotank/derived/gbohner/input/', ...
    'output_folder', '/neurotank/derived/gbohner/output/', ...
    'precomputed_folder', '/neurotank/derived/gbohner/precomputed/', ...
    'data_type', 'frames', ...
    'init_model','filled', ...
    'stabilize', 2, ...
    'niter', 4, ...
    'm', 17, ...
    'mom', 2, ...
    'fig',1, ...
    'spatial_scale',1,...
    'time_scale',1,...
    'mask', 0, ...
    'smooth_filter_mean',2, ...
    'smooth_filter_var',2, ...
    'cells_per_image', 30, ...
    'KS', 4 ...
  );


%For this one no tif images yet, in progress
opt.data_path = '/neurotank/Watkins/2016-02-12/2P/CenterOutReach/site005/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021_Cycle00001_Ch2_000001.ome.tif';

% opt.data_path = '/neurotank/derived/Watkins/2016-02-12/2P_stabilized/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021/image_00001.tif';
% opt.src_string = '*.tif';

%opt.data_path = '/neurotank/derived/Watkins/2016-02-13/2P_stabilized/Tseries_20160213-059/image_00001.tif'; %Stabilizied series
%opt.src_string = '*.tif';

%opt.data_path = '/neurotank/derived/Watkins/2016-02-12/2P_stabilized/Tseries_20160212-018/image_00001.tif';
%opt.src_string = '*.tif';
%opt.timestamp = '20160304T111735';

%opt.data_path = '/neurotank/Watkins/2016-02-13/2P/Exploratory/site002/Tseries_20160213-059/Tseries_20160213-059_Cycle00001_Ch2_000001.ome.tif';
%opt.data_path = '/neurotank/Watkins/2016-02-19/2P/CenterOutReach/site002/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
 %opt.data_path = '/mnt/stanford/neurotank/Watkins/2016-02-19/2P/CenterOutReach/site002/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
 
% % For local run and local data 
 %opt.data_path = '/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
 
 
% opt.input_folder = '~/stanford/input/';
% opt.output_folder = '~/stanford/output/';
% opt.precomputed_folder = '~/stanford/precomputed/';

% Varius param settings
opt.spatial_scale = 0.5;
opt.m = 17;
%opt.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor((opt.m+3)/2-1)); %Should be change when opt.m is changed (%TODO automatically, perhaps with linking to the opt.m variable, symbolic matlab)
%opt.mom = 4;
opt.data_type = 'frames_virtual';
 
[~, opt.file_prefix] = fileparts(opt.data_path);

[opt, ROI_mask, ROIs] = chomp(opt);

get_cell_timeseries;

% load('../data_for_gergo/S1-T52844_masks.mat', 'donut')

% figure(2); imshowpair(ROI_mask, donut>0);