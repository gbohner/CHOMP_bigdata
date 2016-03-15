close all;
clear all;

setenv('CHOMP_ROOT_FOLDER',''); 

cd(fileparts(mfilename('fullpath')));
addpath('./Classes')

%Load an already processed stack (from "input" folder)
load('/neurotank/derived/gbohner/input/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021_Cycle00001_Ch2_000001.ome_20160309T164322.mat');
opt = inp.opt;
opt.init_iter = opt.niter; %If starting with the last iteration, no learning is done, only inference (even if learning is turned on);
opt.learn = 0; %Alternatively just do this and set opt.init_model to 'given' and opt.init_W to already learned basis functions
[opt, ROI_mask, ROIs] = chomp(opt);

timeseries = get_cell_timeseries(opt);