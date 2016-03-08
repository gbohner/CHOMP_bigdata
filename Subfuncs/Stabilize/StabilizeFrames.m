function opt = StabilizeFrames( opt )
%STABILIZEFRAMES Summary of this function goes here
%   Detailed explanation goes here

if ~opt.stabilize,return; end

%For now just implement for frames and frames_virtual
  if strcmp(opt.data_type(1:6),'frames')
    %Align to cascading mean
    if ~exist(get_path(opt, 'raw_stabilized_frames', 1), 'file') || opt.stabilize>1
      [s,mess,messid] = mkdir(get_path(opt, 'raw_stabilized_frames')); %Create folder
      [optimizer, metric]  = imregconfig('monomodal');
      filepath = [fileparts([opt.root_folder opt.data_path]) filesep];
      allfiles = dir([filepath '*' opt.src_string '*']);
      T = size(allfiles,1);
      y0 = imread([filepath allfiles(1).name]);
      y0 = double(zeros(size(y0)));
      %Read the initial images to create an average, always write out the
      %stabilized images
      for i1 = 1:10
        y1 = imread([filepath allfiles(i1).name]);      

        imwrite(y1, get_path(opt, 'raw_stabilized_frames',i1));

        y0 = y0*(i1-1)/i1 + double(y1)./i1;

      end

      %Stabilize the images
      for i1 = 11:T
        y1 = imread([filepath allfiles(i1).name]);     
        y1 = imregister(double(y1),y0,'translation',optimizer,metric);
       
        y0 = y0*(i1-1)/i1 + double(y1)./i1;
        
        imwrite(uint16(y1), get_path(opt, 'raw_stabilized_frames',i1));
        imagesc(y0); colormap gray; axis image; pause(0.1);
      end
    end
    
    opt.data_path = get_path(opt, 'raw_stabilized_frames',1); %Set the data_path to the stabilized frames
    opt.data_path = opt.data_path(length(opt.root_folder)+1:end);
    opt.src_string = '*.tif';
  end

end

