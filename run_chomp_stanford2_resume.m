close all;
clear all;

cd(fileparts(mfilename('fullpath')));
addpath('./Classes')


load('/mnt/stanford/neurotank/derived/gbohner/input/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021_Cycle00001_Ch2_000001.ome_20160307T160319.mat')
opt = inp.opt;
opt.init_iter = opt.niter;
opt.niter = opt.niter + 1;
[opt, ROI_mask, ROIs] = chomp(opt);