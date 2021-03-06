close all;
clear all;

%Runs perfectly on the gatsby leon

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'));

setenv('CHOMP_ROOT_FOLDER','/nfs/data3/gergo/Jim2016/'); %
if ~exist(['./Subfuncs/Compute/Mex/computeGW.' mexext],'file')
  mex('./Subfuncs/Compute/Mex/computeGW.c', '-outdir', './Subfuncs/Compute/Mex/');
end

opt_def_struct = struct(...
    'root_folder', getenv('CHOMP_ROOT_FOLDER'), ...
    'input_folder', 'input/', ...
    'output_folder', 'output/', ...
    'precomputed_folder', 'precomputed/', ...
    'results_folder', 'results/', ...
    'src_string', '*Ch2*', ...
    'data_type', 'frames_virtual', ...
    'init_model',{{'filled','pointlike'}}, ...
    'stabilize', 1, ...
    'niter', 1, ...
    'm', 11, ...
    'mom', 4, ...
    'fig',0, ...
    'spatial_scale',0.5,...
    'time_scale',1,...
    'mask', 0, ...
    'cells_per_image', 50, ...
    'KS', 10 ...
  );


%Initialize all datasets to the same options (shouldn't use the same
%chomp_options object and deal(obj) because it is passed by reference then)
opts = cell(5,1);
for n = 1:numel(opts)
    opts{n} = chomp_options(opt_def_struct);
end

cellsize  = {};
 
%Set the individually differing parameters (possibly automatically later)
opts{1}.data_path = 'vStim-043/vStim-043_Cycle00001_CurrentSettings_Ch2_000001.tif';
[~, opts{1}.file_prefix] = fileparts(opts{1}.data_path);
%Make sure apparent neuron sizes are the same, as well as they are within
%the used basis function size
%cell size 25

opts{2}.data_path = 'vStim-044/vStim-044_Cycle00001_CurrentSettings_Ch2_000001.tif';
[~, opts{2}.file_prefix] = fileparts(opts{2}.data_path);
%Make sure apparent neuron sizes are the same, as well as they are within
%the used basis function size
%cell size 25


opts{3}.data_path = '/visualStim-002/visualStim-002_Cycle00001_CurrentSettings_Ch1_000001.tif';
[~, opts{3}.file_prefix] = fileparts(opts{3}.data_path);
%Make sure apparent neuron sizes are the same, as well as they are within
%the used basis function siz
%cell size 25

opts{4}.data_path = 'visualStim-004/visualStim-004_Cycle00001_CurrentSettings_Ch1_000001.tif';
[~, opts{4}.file_prefix] = fileparts(opts{4}.data_path);
%Make sure apparent neuron sizes are the same, as well as they are within
%the used basis function siz
%cell size 25

opts{5}.data_path = '/visualStim-005/visualStim-005_Cycle00001_CurrentSettings_Ch1_000001.tif';
[~, opts{5}.file_prefix] = fileparts(opts{5}.data_path);
%Make sure apparent neuron sizes are the same, as well as they are within
%the used basis function size
%cell size 25

%%

%Initialize basis subspaces
W_cur = Model_initialize(opts{1});

gtic = tic;

maxiter = 11;
for iters = 1:maxiter
  %Run chomp one inference step forward for each dataset
  fprintf('Running batch inference on %d datasets, starting iteration %d/%d...\n',numel(opts),iters,maxiter);
  parfor n = 1:numel(opts)
    opts{n}.niter = iters;
    opts{n}.init_iter = iters-1;
    opts{n}.init_W = W_cur;
    %For last iteration do inference on many proposed objects
    if iters==maxiter
      opts{n}.cells_per_image = opts{n}.cells_per_image * 10;
      opts{n}.learn = 0;
    end
    opts{n} = chomp(opts{n});
  end
  
  %For last iteration no learning
  if iters == maxiter, break; end;
  
  fprintf('Running batch learning on %d datasets, starting iteration %d/%d...\n',numel(opts),iters,maxiter);
  %Do the learning of new W jointly for all datasets
  for n=1:numel(opts)
    inp = load(get_path(opts{n}));
    datas{n} = inp.inp.data;
    outp = load(get_path(opts{n},'output_iter',opts{n}.niter),'model');
    Hs{n} = outp.model.H;
  end
  for type = 1:opts{1}.NSS, Wblocked{type} = W_cur(:,opts{1}.Wblocks{type}); end
  parfor type = 1:opts{1}.NSS
    Wblocked{type} = update_dict(datas,Hs,Wblocked{type},opts,iters+1,type);
  end
  for type = 1:opts{1}.NSS, W_cur(:,opts{1}.Wblocks{type}) = Wblocked{type}; end
  
  fprintf('Finished batch iteration %d/%d, time elapsed is %.2f\n',iters,maxiter,toc(gtic));
end 

%%

%Get time series for all datasets
timeseries = cell(numel(opts),1);
parfor n = 1:numel(opts)
  timeseries{n} = get_cell_timeseries(opts{n});
end



