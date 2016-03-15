close all;
clear all;

%Example local run

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'))

setenv('CHOMP_ROOT_FOLDER', '');

opt = chomp_options(); %Initialize default options

% % Setup folder structure for local run and local data 
opt.data_path = '~/stanford/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001/Tseries_20160219_Watkins_CenterOutReach_time20160219.133302.428-001_Cycle00001_Ch2_000001.ome.tif';
opt.src_string = '*Ch2*';

opt.input_folder = '~/stanford/input/';
opt.output_folder = '~/stanford/output/';
opt.precomputed_folder = '~/stanford/precomputed/';
opt.results_folder = '~/stanford/results/';
 

[~, opt.file_prefix] = fileparts(opt.data_path); % Important if you want to have nice file prefixes (corresponding to folder name)
[opt, ROI_mask, ROIs] = chomp(opt);
get_cell_timeseries(opt);