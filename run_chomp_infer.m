close all;
clear all;

setenv('CHOMP_ROOT_FOLDER',''); 

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'));

%Load an already learned model (from "output" folder, watch out for correct timestamp)
load('~/stanford/output/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome_20160315T004143_iter_4.mat');

%Use the same opt struct to ensure similar preprocessing, but change data
%path, spatial_scale and remove the timestamp
opt = chomp_options(model.opt.export_struct); %New constructor to make sure we don't overwrite the original instance (even though here it would work ok, it matters for batch processing)
opt.data_path = '~/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
[~, opt.file_prefix] = fileparts(opt.data_path); %Set the new file prefix
opt.spatial_scale = 1; %Set it to match the same expected cell size as in the previous dataset
opt.timestamp = '';

%Give the learned basis functions as input (specify we gave it for each
%object type)
for type = 1:opt.NSS, opt.init_model{type} = 'given'; end
opt.init_W = model.W; % It is possible to only give a few, and reinitialize other cell types, but for inference we give all

%Set up the iterations for learning
opt.init_iter = 0;
opt.niter = 1;
opt.learn = 0;

%Set the number of objects CHOMP will return, then you can cut the non-important ones later
opt.cells_per_image = 200; 

%Run the inference
[opt, ROI_mask, ROIs] = chomp(opt);

%Get time series
timeseries = get_cell_timeseries(opt);