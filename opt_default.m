% Set up options
function opt = opt_default()



opt = struct(); 

opt.code_path = [fileparts(mfilename('fullpath')) filesep]; %Package directory
opt.data_path = 'default_path'; %Input data file (stack, or initial frame)
opt.input_folder = './tmp/input/'; %Store preprocessed input data, if you change it, use full path
opt.output_folder = './tmp/output/'; %Output folder, if you change it, use full path
opt.precomputed_folder = './tmp/precomputed/'; % Stores precomputed tensors
opt.file_prefix = 'test'; %Prefix for file names

% Model setup
opt.m = 13; % Basis function size in pixels
opt.NSS = 1; % Number of object types
opt.KS = 4; % Dimensionality of space per object type (i.e. number of basis functions per object type)
opt.init_model = 'donut'; % 'filled', 'donut', 'pointlike', \\ %TODO: 'supervised', 'multi'
opt.init_W = [];

% Data extraction and preprocessing
opt.spatial_scale = 1; % Rescale data spatially (so that cell size matches basis function size)
opt.time_scale = 0.15; % Rescale data temporally
opt.whiten = 1;
opt.smooth_filter_mean = opt.m; %smoothing filter size for mean image
opt.smooth_filter_var = opt.m; %smoothing filter size for variance
opt.data_type = 'frames'; %Input data type (frames / stack / json / matxyt)
opt.src_string = 'Ch2_*'; %in case of loading multiple frames from a directory, look for this substring to load files (choose channel eg)
opt.mask = 0; % Set if the region of interest is only part of the image stack.
opt.mask_overwrite = 0; %If a previous intermediate file already has mask, do you want the to create a new one (1) or automatically use the old one (0).
opt.mask_image = []; % you can input your own binary mask image if needed
opt.mom = 1; %Number of moments used

% Learning parameters
opt.niter = 3; % number of iterations
%opt.relweight = 10; % weighting between importance of covariance / mean (automatically set to 'optimal' value in Shared_main/extract_coefs.m)
opt.fig = 1; %Whether to visualize or not during learning
opt.ex          = 1; % what example image to display during training
opt.cells_per_image = 100; % a rough estimate of the average number of cells per image
opt.relweight = 1; %Relative weight between mean and correlation coeff.
opt.MP      = 0; % somewhat redundant: if set to 1 always uses one subspace per object
opt.inc     = 2; % every opt.inc iterations estimate a new subspace
opt.warmup = 1;
opt.learn   = 1; % do learning?
opt.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor(opt.m/2-1)); % Specified distance based function (leave as [] if not desired)
opt.learn_decomp = 'HOSVD'; % HOSVD or MTF (MTF not implemented yet, %TODO - write R wrapper to use Kahn2015 code)


% Extracting ROIs
opt.ROI_type = 'quantile';
opt.ROI_params = [0.7];


% Misc parameters
opt.cleanup = 0;

opt.timestamp = '20151112T173059';


end