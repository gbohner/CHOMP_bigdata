close all;
clear;

%Example local run

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'))

setenv('CHOMP_ROOT_FOLDER', '');

opt = chomp_options(); %Initialize default options

% % Setup folder structure for local run and local data 
%opt.data_path = '~/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
opt.data_path = '/Users/gergobohner/Dropbox/Gatsby/Research/forUCL/stacksForUCL/GreenChan/GreenChanB_0000.tif';
%opt.data_path = '/Users/gergobohner/Dropbox/Gatsby/Research/forUCL/stacksForUCL/RedChan/RedChanD_0000.tif';

opt.src_string = '*.tif';
opt.stabilize = 0;
opt.spatial_scale = 1;
opt.NSS = 1;
opt.KS = 8;
opt.init_model = 'filled';
opt.niter = 8;

% opt.input_folder = '~/stanford/input/';
% opt.output_folder = '~/stanford/output/';
% opt.precomputed_folder = '~/stanford/precomputed/';
% opt.results_folder = '~/stanford/results/';
 
opt.fig = 1;

[~, opt.file_prefix] = fileparts(opt.data_path); % Important if you want to have nice file prefixes (corresponding to folder name)

%%
%Learning with 100 cells;
opt.cells_per_image = 60;
[opt, ROI_mask, ROIs] = chomp(opt);

%%
%Inferring with 500 cells;
opt.cells_per_image = 500;
opt.init_iter = opt.niter; %Just running the inference
opt.niter = opt.niter+1; %Just running the inference
opt.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor(opt.m/2-1)); %@(grid_dist, sharp)logsig(sharp*grid_dist-floor(sharp*2*obj.m/2-1));
[opt, ROI_mask, ROIs] = chomp(opt);


get_cell_timeseries(opt);