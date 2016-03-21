classdef chomp_options < handle
  %OPTIONS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetObservable, GetObservable, AbortSet)
    root_folder = ''; % Can set a root folder in case of remote work (sshfs or runs from different users/systems mac vs linux directory structure). All other folders below will be calculated relative from the root folder
  end
  
  properties
    %Setup folder structure
     code_path = [fileparts(mfilename('fullpath')) filesep]; %Package directory
     data_path = 'default_path'; %Input data file (stack, or initial frame)
     input_folder = './tmp/input/'; %Store preprocessed input data, if you change it, use full path
     output_folder = './tmp/output/'; %Output folder, if you change it, use full path
     precomputed_folder = './tmp/precomputed/'; % Stores precomputed tensors
     results_folder = './tmp/results/'; %Store the extracted ROIs and timeseries
     file_prefix = 'test'; %Prefix for file names
 
     % Model setup
     m = 17; % Basis function size in pixels
     NSS = 2; % Number of object types
     KS = 4; % Dimensionality of space per object type (i.e. number of basis functions per object type)
     init_model = {'filled', 'pointlike'}; % 'filled', 'donut', 'pointlike', \\ %TODO: 'supervised', 'multi'
     init_W = [];

     % Data extraction and preprocessing
     stabilize = 1;
     spatial_scale = 1; % Rescale data spatially (so that cell size matches basis function size)
     time_scale = 1; % Rescale data temporally
     whiten = 1;
     smooth_filter_mean % = m; %smoothing filter size for mean image
     smooth_filter_var % = m; %smoothing filter size for variance
     data_type = 'frames_virtual'; %Input data type (frames / stack / json / matxyt)
     src_string = 'Ch2_*'; %in case of loading multiple frames from a directory, look for this substring to load files (choose channel eg)
     mask = 0; % Set if the region of interest is only part of the image stack.
     mask_overwrite = 0; %If a previous intermediate file already has mask, do you want the to create a new one (1) or automatically use the old one (0).
     mask_image = []; % you can input your own binary mask image if needed
     mom = 2; %Number of moments used
     A % Mean image for spatial whitening
     B % Variance for spatial whitening

     % Learning parameters
     niter = 4; % number of iterations
     init_iter = 0; %Set the initial iteration. If not 0, search for the file with appropriate name
      %     relweight = 10; % weighting between importance of covariance / mean (automatically set to 'optimal' value in Shared_main/extract_coefs.m)
     fig = 0; %Whether to visualize or not during learning
     cells_per_image = 30; % the maximum number of objects to infer
     warmup = 1;
     learn   = 1; % do learning?
     spatial_push % = @(grid_dist)logsig(0.5*grid_dist-floor(options.m/2-1)); % Specified distance based function (leave as [] if not desired)
     learn_decomp = 'COV'; % COV, LMSVD, HOSVD or MTF (MTF not implemented yet, %TODO - write R wrapper to use Kahn2015 code)


    % Extracting ROIs
     ROI_type = 'quantile_dynamic_origsize';
     ROI_params = [0.6];


     % Misc parameters
     
     cleanup = 0;

     timestamp % Timestamping
     
     
  end
  
  properties (Dependent)
    Wblocks
  end
  
  methods
    function obj = chomp_options(varargin)
      %Object constructor
      
      %Construct from a struct
      if nargin==1
        if isa(varargin{1},'struct')
          fns = fieldnames(varargin{1});
          for i1 = 1:numel(fns)
            try
              obj.(fns{i1}) = varargin{1}.(fns{i1});
            catch ME
              rethrow(ME); %TODO better error message
            end
          end
        else
          error('CHOMP:chomp_options:bad_constructor','Wrong class constructor call for chomp_options');
        end
      end
      
          
      %Construct from ('field', 'value') argument pairs
      i1 = 1;
      while i1<nargin
        obj.(varargin{i1}) = varargin{i1+1};
        i1 = i1+2;
      end
      
      %Compute the derived properties
      obj = obj.derive_from_m(); 
      
      obj = obj.assert();
    end
    
    function obj = assert(obj)
      %Check for specific parameters to be in the correct format;
      if ~iscell(obj.init_model), obj.init_model = {obj.init_model}; end
      assert(numel(obj.init_model)==obj.NSS, 'CHOMP: Object type # discrepency');
    end
    
    
    function obj = derive_from_m(obj)
      %On construct, makes sure that certain properties are set correctly
      %in relation to basis function (i.e. expected cell) size
      obj.smooth_filter_mean = obj.m;
      obj.smooth_filter_var = obj.m;
      obj.spatial_push = @(grid_dist)logsig(0.5*grid_dist-floor(obj.m/2-1)); %@(grid_dist, sharp)logsig(sharp*grid_dist-floor(sharp*2*obj.m/2-1));
    end
    
    function s = export_struct(obj, varargin)
      p = properties(obj);
      for i1 = 1:numel(p)
        s.(p{i1}) = obj.(p{i1});
      end
    end

    function blocks = get.Wblocks(obj)
      blocks = cell(obj.NSS, 1);
      for type = 1:obj.NSS
        blocks{type} = ((type-1)*obj.KS+1):(type*obj.KS);
      end
    end
    
    function set.Wblocks(obj, val)
      %Just to suppress errors coming from no set method, it doesn't do
      %anything.
    end
  end
  
  
  
  
end

