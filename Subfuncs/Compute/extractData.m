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
%       user_ans = input('Continue processing with an updated copy of the intermediate file (y) or load the old one (n)? y/n  ', 's');
%       if strcmp('y', user_ans) % get new time stamp, copy the file, update options
%         inp.opt = struct_merge(inp.opt, opt);
%         intermediate_path = get_path(inp.opt);
%         save(intermediate_path, 'inp', '-append'); %overwrite the old options file within the file with the new time stamp
%       end
    end
% TODO: fix the mask update part    
%     if opt.mask %Check if we need to update the user mask in the new file
%       if ~isempty(opt.mask_image) 
%         UserMask = opt.mask_image;
%         save(intermediate_path, 'UserMask', '-append');
%       elseif opt.mask_overwrite==1 || inp.opt.mask==0
%         load(intermediate_path,'y');
%         figure(1); imagesc(y); axis square; colormap gray;
%         UserMask = roipoly(mat2gray(y));
%         save(intermediate_path, 'UserMask', '-append');
%       end
%     end
  end
     
  
else %Handle the data loading-preprocessing-saving

  %Set the input data path
  data_path = [opt.root_folder opt.data_path];
  
  if ~exist(data_path,'file')
      error('CHOMP:preprocess:noinputdata', ...
      ['There is no file at the given input destination to read']);
  end
  
  %Read the data into variable Is (image-stack)
  if strcmp(opt.data_type, 'stack')
    info = imfinfo(data_path);
    T = numel(info); % Number of frames
    sz = size(imresize(imread(data_path,1),opt.spatial_scale)); % Image size
    if T>1
      data.proc_stack.Y = zeros([sz(1), sz(2), T]); % Image stack
      for i2 = 1:T
          data.proc_stack.Y(:,:, i2) = imresize(double(imread(data_path, i2)),opt.spatial_scale);
      end
    else
      data.proc_stack.Y = imresize(double(imread(data_path)),opt.spatial_scale);
    end
    %Store the file link to the image stack
    data.raw = datastore(data_path,'Type', 'image');
  elseif strcmp(opt.data_type, 'frames')
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
    
  elseif strcmp(opt.data_type, 'json')
      conf = loadjson(data_path); %For json stuff give the configuration file as input path
      filepath = [fileparts(data_path) filesep];
      allfiles = dir([filepath '*' opt.src_string '*']);
      fid = fopen([filepath allfiles(1).name],'r');
      
      sz = size(imresize(double(reshape(fread(fid, 'uint16'), conf.dims(1), conf.dims(2))),opt.spatial_scale,'bicubic'));
      fclose(fid);
      T = size(allfiles,1);
      data.proc_stack.Y = zeros([sz(1), sz(2), T]);
      for i2 = 1:T
        fid = fopen([filepath allfiles(i2).name],'r');
        data.proc_stack.Y(:,:,i2) = imresize(double(reshape(fread(fid, 'uint16'), conf.dims(1), conf.dims(2))),opt.spatial_scale,'bicubic');
        fclose(fid);
      end
      %TODO data.raw
  elseif strcmp(opt.data_type, 'matxyt')
      Ytmp = load(data_path);
      Ytmp = double(Ytmp.data);
      sz = size(imresize(Ytmp(:,:,1), opt.spatial_scale, 'bicubic'));
      T = size(Ytmp,3);
      data.proc_stack.Y = zeros([sz(1), sz(2), T]);
      for i2 = 1:T
        data.proc_stack.Y(:,:,i2) = imresize(Ytmp(:,:,i2), opt.spatial_scale, 'bicubic');
      end
      data.raw = matfile(data_path);
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
      s2 = 0;
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
          disp(s2)
        end
      end
  end % data reading
  
  %Downsample in time (average in subsequent time windows)
  szProc = chomp_size(data.proc_stack,'Y');
  tmp = zeros([szProc(1:2),floor(T *opt.time_scale)]); %TODO make as temporary file maybe
  if opt.time_scale < 1
    T = floor(T *opt.time_scale);
    parfor i2 = 1:T
      tmp(:,:,i2) = mean(data.proc_stack.Y(:,:,ceil((i2-1)/opt.time_scale)+1:floor(i2/opt.time_scale)),3); %TODO - do this parallel?
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
  tmp = zeros(chomp_size(data.proc_stack,'Y')); %TODO
  if opt.whiten
      parfor t1 = 1:T
        tmp(:,:,t1) = squeeze(data.proc_stack.Y(:,:,t1)-opt.A) ./ (opt.B.^0.5);
      end
  end
  data.proc_stack.Y = chomp_data(data.proc_stack.Y.Source,tmp);

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
  
  
  

