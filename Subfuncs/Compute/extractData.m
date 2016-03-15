function opt = extractData( opt )
%EXTRACTDATAFROMTIF Extracts the relavant information from an input tif
%stack 
% - give input tif file name
% - output path for large .mat file
% - other arguments:
%   - options struct
%   -  ...

  opt=opt; %avoids matlab error of "output not assigned..."

  %Set the path for the intermediate preprocessed "input" file
  intermediate_path = get_path(opt);
  
%Check if we want to just load an already preprocessed file
if exist(intermediate_path, 'file')
  %If it already exists, check if the preprocessing options are the same, and just load the files
  inp = load(intermediate_path, 'inp');
  inp = inp.inp;
  
  if ~struct_contain(opt, inp.opt, {'spatial_scale', 'time_scale', 'whiten', 'smooth_filter_mean' , 'smooth_filter_var'})
    error('CHOMP:preprocess:outdatedoptions', ...
      ['The new option struct has different preprocessing options than the', ...
       'intermediate file you want to use it with, consider removing manual timestamp' ...
       'from your input options file to create a new preprocessed file']);
  else
    if ~struct_contain(opt, inp.opt) %some settings has changed, give warning, and let user create a new copy of the file with new timestamp
      warning('CHOMP:preprocess:outdatedoptions', 'The new option struct has slightly different options than the intermediate file you want to use it with.');
      inp.opt = struct_merge(inp.opt, opt);
      intermediate_path = get_path(inp.opt);
      save(intermediate_path, 'inp', '-append'); %overwrite the old input file
    end
  end
     
  
else %Handle the data loading-preprocessing-saving

  %Set the input data path
  data_path = [opt.root_folder opt.data_path];
  
  if ~exist(data_path,'file')
      error('CHOMP:preprocess:noinputdata', ...
      ['There is no file at the given input destination to read']);
  end
  
  %Read the data into a preprocessed and a raw stack, as well as store the
  %original files in a datastore
  if strcmp(opt.data_type, 'frames')
    filepath = [fileparts(data_path) filesep];
    allfiles = dir([filepath '*' opt.src_string '*']);
    T = size(allfiles,1);
    sz = size(imresize(imread([filepath allfiles(1).name]),opt.spatial_scale));
    %Store the file links to the raw data
    data.raw = datastore(strcat(filepath, {allfiles.name}),'FileExtension','.tif','Type', 'image');
    raw_stack_done = 0; if exist(get_path(opt,'raw_virtual_stack'),'file'), raw_stack_done = 1; end;
    if raw_stack_done, data.raw_stack.Y = chomp_data(get_path(opt,'raw_virtual_stack')); 
    else
      data.raw_stack.Y = chomp_data(get_path(opt,'raw_virtual_stack'), double(imread([filepath allfiles(1).name])));
    end
    %Save the preprocessed stack
    data.proc_stack.Y = zeros([sz(1), sz(2), T]); % Image stack
    im_cur = double(imread([filepath allfiles(1).name]));
    data.proc_stack.Y(:, :,1) = imresize(im_cur,opt.spatial_scale,'bicubic');
    for i2 = 2:T
      im_cur = double(imread([filepath allfiles(i2).name]));
      data.proc_stack.Y(:, :,i2) = imresize(im_cur,opt.spatial_scale,'bicubic');
      if ~raw_stack_done, append(data.raw_stack.Y, im_cur); end
    end
  elseif strcmp(opt.data_type, 'frames_virtual')
      filepath = [fileparts(data_path) filesep];
      allfiles = dir([filepath '*' opt.src_string '*']);
      T = size(allfiles,1);
      sz = size(imresize(imread([filepath allfiles(1).name]),opt.spatial_scale));
      %Store the file links to the raw data
      data.raw = datastore(strcat(filepath, {allfiles.name}),'FileExtension','.tif','Type', 'image');
      raw_stack_done = 0; if exist(get_path(opt,'raw_virtual_stack'),'file'), raw_stack_done = 1; end;
      if raw_stack_done, data.raw_stack.Y = chomp_data(get_path(opt,'raw_virtual_stack')); 
      else
        data.raw_stack.Y = chomp_data(get_path(opt,'raw_virtual_stack'), double(imread([filepath allfiles(1).name])));
      end
      
      Ytmp = zeros([sz(1:2),100]);
      %Minibatch read and dump to the matfile variable
      s2 = 0; charcount = 0;
      for i2 = 1:T
        im_cur = double(imread([filepath allfiles(i2).name]));
        if (~raw_stack_done) && (i2>1), append(data.raw_stack.Y, im_cur); end
        Ytmp(:,:,mod(i2,100)+floor(i2/100)*100) = imresize(im_cur,opt.spatial_scale,'bicubic');
        if (mod(i2,100) == 0) || (i2 == T)
          if s2==0
            data.proc_stack.Y = chomp_data(get_path(opt,'virtual_stack'),Ytmp(:,:,1:(mod(min(((s2+1)*100),T)-1,100)+1)));
          else
            append(data.proc_stack.Y, Ytmp(:,:,1:(mod(min(((s2+1)*100),T)-1,100)+1)));
          end
          s2 = s2+1;
          if charcount>0, for c1 = 1:charcount, fprintf('\b'); end; end
          charcount = fprintf('Reading images and creating virutal stack... %d/%d',s2*100, T);
        end
      end
  end % data reading
  
  
  %Downsample in time (average in subsequent time windows)
  szProc = chomp_size(data.proc_stack,'Y');
  if opt.time_scale < 1
    T = floor(T *opt.time_scale);
    for t1 = 1:T
      tmp = mean(data.proc_stack.Y(:,:,ceil((t1-1)/opt.time_scale)+1:floor(t1/opt.time_scale)),3);
      %Store results in place
      if strfind(opt.data_type, 'virtual')
        overwrite_frame(data.proc_stack.Y, tmp, t1);
      else
        data.proc_stack.Y(:,:,t1) = tmp;
      end  
    end
    data.proc_stack.Y = chomp_data(data.proc_stack.Y.Source,tmp);
  end
  
  %Get mean image of the raw data
  y_orig = double(readimage(data.raw,1))./numel(data.raw.Files);
  for i1 = 2:numel(data.raw.Files)
    y_orig = y_orig + double(readimage(data.raw,i1))./numel(data.raw.Files);
  end
 
  %Get mean image and variance, for cycle for virtual stack handling
  y = data.proc_stack.Y(:,:,1)./T;
  for i1 = 2:T
    y = y + squeeze(data.proc_stack.Y(:,:,i1))./T;
  end
  
  %Get the variance
  if T>1
    V = zeros(size(y));
    for i1 = 1:T
      V = V + (squeeze(data.proc_stack.Y(:,:,i1))-y).^2./T;
    end
  else
    V = ones(size(y));
  end
  
%   %Fix edge_effect problems later %TODO
%   edge_effect = conv2(ones(size(y)),ones(opt.m),'same');
  
  % Apply normalizing filters to mean image for visualization purposes and
  % to get "mean spatial mean (A) and mean spatial variance (B)"
  if T>1
    [y, A, B] = normal_img(double(y), opt.smooth_filter_mean , opt.smooth_filter_var ,V);
  else
    [y, A, B] = normal_img(double(y), opt.smooth_filter_mean , opt.smooth_filter_var);
  end
  opt.A = A;
  opt.B = B;
  
  %Apply whitening to the whole stack (whitening with respect to the
  %smoothing filter sizes set, not on a pixel-by-pixel level!)
  
  fprintf('\nPreprocessing image stack...');
  if opt.whiten
      for t1 = 1:T
        tmp = squeeze(data.proc_stack.Y(:,:,t1)-opt.A) ./ (opt.B.^0.5);
        if strfind(opt.data_type, 'virtual')
          overwrite_frame(data.proc_stack.Y, tmp, t1);
        else
          data.proc_stack.Y(:,:,t1) = tmp;
        end       
      end
  end
  

  inp = chomp_input(opt,data, y, y_orig,V); %The input class that stores raw and preproc data
  
  
   %let the user create a mask over the image to limit the processed area
  if opt.mask
    if ~isempty(opt.mask_image) 
        inp.UserMask = opt.mask_image;
    else
      imshow(y);
      inp.UserMask = roipoly(mat2gray(y));
    end
  end
  
  save(intermediate_path, 'inp','-v7.3');
  
  

end

 
  
  clearvars -except 'opt'

end
  
  
  

