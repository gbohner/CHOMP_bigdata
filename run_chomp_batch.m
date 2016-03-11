close all;
clear all;

cd(fileparts(mfilename('fullpath')));
addpath('./Classes')

setenv('CHOMP_ROOT_FOLDER',''); %

opt_def = chomp_options(...
    'root_folder', getenv('CHOMP_ROOT_FOLDER'), ...
    'input_folder', '/neurotank/derived/gbohner/input/', ...
    'output_folder', '/neurotank/derived/gbohner/output/', ...
    'precomputed_folder', '/neurotank/derived/gbohner/precomputed/', ...
    'data_type', 'frames_virtual', ...
    'init_model','filled', ...
    'stabilize', 1, ...
    'niter', 4, ...
    'm', 17, ...
    'mom', 2, ...
    'fig',1, ...
    'spatial_scale',1,...
    'time_scale',1,...
    'mask', 0, ...
    'cells_per_image', 15, ...
    'KS', 6 ...
  );


%Initialize all datasets to the same options
opts = cell(2,1);
[opts{:}] = deal(opt_def);


%Set the individually differing parameters (possibly automatically later)
opts{1}.data_path = '/neurotank/Watkins/2016-02-12/2P/CenterOutReach/site005/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-021_Cycle00001_Ch2_000001.ome.tif';
%Make sure apparent neuron sizes are the same, as well as they are within
opts{1}.spatial_scale = 0.7;

opts{2}.data_path = '/neurotank/Watkins/2016-02-12/2P/CenterOutReach/site005/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-024/Tseries_20160212_Watkins_CenterOutReach_time20160212.123454.112-024_Cycle00001_Ch2_000001.ome.tif';
opts{2}.spatial_scale = 1.4;

W_cur = Model_initialize(opts{1});

for iters = 1:5
  %Run chomp one inference step forward for each dataset
  parfor n = 1:numel(opts)
    opts{n}.niter = iters;
    opts{n}.init_iter = iters-1;
    opts{n}.W_init = W_cur;
    opts{n} = chomp(opts{n});
  end
  
  %Do the learning of new W jointly for all datasets
  parfor n=1:numel(opts)
    inp = load(get_path(opts{n}))
    datas{n} = inp.inp.data;
    outp = load(get_path(opts{n},'output_iter',opts{n}.niter),'model')
    Hs{n} = outp.model.H;
  end
  W_cur = update_dict(datas,Hs,W_cur,opts,iters+2);
end



