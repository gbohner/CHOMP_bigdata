function out = readStack(path, varargin)
%READSTACK Reads a chunk of the data stack into an output array
%path to file, and optionally: frames, patches
  p = inputParser();
  p.addRequired('path',@(x)exist(x,'file'))
  p.parse(path);
  
  fid = fopen(p.Results.path,'r');
  dims = fread(fid, 1, 'double');
  szData = fread(fid,uint16(dims),'double')';
  headerStr = char(fread(fid,30,'char'))';
  
  p.addOptional('frames',1:szData(end),@isfloat);
  p.addParameter('patch',[],@(x)any([isempty(x),isstruct(x),iscell(x)]));
  p.addParameter('minDiskAccess',0,@isnumeric);
  p.parse(path, varargin{:})
  
  %disp(p.Results)
  
  
  %Determine the shape of the output array (that is frame or patchsize x num_frames x num_patches)
  if isempty(p.Results.patch)
    %Return all the frames required
    szOut = [szData(1:end-1), numel(p.Results.frames), 1];
  elseif isstruct(p.Results.patch)
    szOut = [numel(p.Results.patch.x), numel(p.Results.patch.y), numel(p.Results.frames), 1];
  elseif iscell(p.Results.patch)
    %TODO check all patches and return the maximum one, or look for out-of
    %-dataset indexing within this file
    szOut = [numel(p.Results.patch{1}.x), numel(p.Results.patch{1}.y), numel(p.Results.frames), numel(p.Results.patch)];
  end
  
  %TODO Check if szOut amount of data fits into memory
  
  %preallocate the output array
  %szOutProd = [prod(szOut(1:2)),szOut(3:4)];
  out = zeros(szOut); %simpler indexing
  
  %Determine the set of bytes we want to read, and where to put them within
  %the output array (assignment tensor)
  if p.Results.minDiskAccess && iscell(p.Results.patch)
    %PPrecompute the absolute minimal set of bytes we need to read
    %First check for a single frame which bytes to read and where to write
    %them
    frameBytes = zeros([szData(1:end-1),1]); %Which ones to read
    ByteWriteTo = sparse(prod(szData(1:2)),prod([szOut(1:2),szOut(4)])); %Where to write assignment tensor %TODO
    for c1 = 1:numel(p.Results.patch)
      frameBytes(p.Results.patch{c1}.x,p.Results.patch{c1}.y) = 1;
    end

    frameBytes = frameBytes(:);

    frameByteSkips = 0:(numel(frameBytes)-1); %How much to skip before each read
    frameByteSkips(frameBytes==0) = [];
    frameByteSkips = [frameByteSkips(1) diff([frameByteSkips,numel(frameBytes)])-1]; %Include skipping to the very last one
  end
  
  frameSkips = [0, p.Results.frames]; %How many frames to skip before each read one
  frameSkips = diff(frameSkips)-1;
  %Do the reading procedure frame by frame, fill up the output array
  %accordingly
  
  numSkipAmount = 8; % 8 bytes per double
  frameSkipAmount = numSkipAmount * prod(szData(1:end-1)); %single frame
  
  
  if isstruct(p.Results.patch)
    rowSkipBef = (min(p.Results.patch.x)-1)*numSkipAmount;
    rowRead = numel(p.Results.patch.x);
    rowSkipAft = szData(1)-max(p.Results.patch.x)*numSkipAmount;
    colSkipBef = (min(p.Results.patch.y)-1)*szData(1)*numSkipAmount;
    colSkipAft = (szData(2)-max(p.Results.patch.y))*szData(1)*numSkipAmount;
  end
  
  all_patches = 1:szOut(end);
  
  for t1 = 1:numel(p.Results.frames)
    %Skip frames
    fseek(fid,frameSkips(t1)*frameSkipAmount,'cof');
    if isempty(p.Results.patch)
      out(:,:,t1) = reshape(fread(fid,prod(szOut(1:2)),'double'),szOut(1:2));
    elseif isstruct(p.Results.patch) %TODO still multiple possibilities, if frameBytes is dense, maybe worth loading the full frame
      fseek(fid,colSkipBef,'cof');
      for j1 = 1:szOut(2) %Read columns of the patch
        fseek(fid,rowSkipBef,'cof');
        out(:,j1,t1,1) =  fread(fid,rowRead,'double');
        fseek(fid,rowSkipAft,'cof');
      end
      fseek(fid,colSkipAft,'cof'); %Arrive at start of next frame
    elseif iscell(p.Results.patch)
      if ~p.Results.minDiskAccess %Read full frames and extract the required patch blocks
      elseif p.Results.minDiskAccess %absolutely minimal disk access, sacrificing a bit more computational speed via indexing
        %Skip bytes
        numCount = 0;
        withinPatchCount = zeros(szOut(end),1); %TODO: if some rows or whatever is truncated, adjust accordingly maybe with a withinPatchCountSkips
        for i1 = 1:(numel(frameByteSkips)-1)
          fseek(fid,frameByteSkips(i1)*numSkipAmount,'cof');
          numCount = numCount+frameByteSkips(i1)+1;
          patchWritten = d2b(frameBytes(numCount), szOut(end)); %This determines which patches I am writing the number into
          withinPatchCount = withinPatchCount + patchWritten;
          cur_patches = all_patches(patchWritten);
          cur_val = fread(fid,1,'double');
          for j1 = cur_patches
            out(withinPatchCount(j1),t1,j1) = cur_val;
          end
        end
        fseek(fid,frameByteSkips(end)*numSkipAmount,'cof');  
      end
    end
  end

  fclose(fid);

  out = reshape(out,szOut);

end

