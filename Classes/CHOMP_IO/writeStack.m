function writeStack( path, data, varargin )
      %Writes the data into binary file(s)
      %Setup: First 'double' number is ndims of data
      % Then the next ndims 'double' numbers is individual dimensions of data
      % Then the next 100 'char' header is:
      %   1-10: The number format used
      %  11-80: Arbitrary metadata (file prefix)
      % 86-100: The timestamp of the dataset
      % From there on we store the individual frames in order:
      % Columns x Rows x Frames (equivalent to Matlab fwrite(fid,
      % data(:),number_format) format);
      
     
      
      %File format:
      p = inputParser();
      p.addRequired('path',@ischar)
      p.addRequired('data',@isnumeric)
      p.addParameter('prefix','DEFAULT_PREFIX',@ischar)
      p.addParameter('timestamp',datestr(now, 30), @ischar)
      p.addParameter('append',0,@(x)all([isnumeric(x),exist(path,'file')])); %Appends to the end of the file
      p.addParameter('overwrite_frame',0,@(x)all([isnumeric(x),exist(path,'file')])); %overwrite frame #overwrite_frame
      p.addParameter('truncateToFrame',0,@(x)all([isnumeric(x),exist(path,'file')])); %
      p.addParameter('number_format','double',@ischar);
      p.parse(path, data,varargin{:})
      
      %Get the number format we are going to use
      number_format = p.Results.number_format;
      tmp = 1; tmp = cast(tmp,number_format); tmp = whos('tmp');
      number_format_bytes = tmp.bytes;
      
      to_modify = any([p.Results.append, p.Results.overwrite_frame, p.Results.truncateToFrame]);
      
      %Open the file at the given location
      if to_modify
        fid = fopen(p.Results.path,'rb+');     
      else
        fid = fopen(p.Results.path,'w');
      end
      
      %Convert the data to required number format
      data = cast(p.Results.data,number_format);
      
      %Write into the header [ndims(data), size(data)]
      if to_modify
        %Get the size of already existing data
        frewind(fid);
        dims = fread(fid, 1, 'double');
        szData = fread(fid,uint16(dims),'double')';
        orig_number_format = strtrim(char(fread(fid,10,'char')'));
        %Assert that dimensions and the number_format is correct
        if ~p.Results.truncateToFrame, assert(all(szData(1:2)==[size(data,1),size(data,2)]),'Dimension for appending are inconsistant'); end
        assert(strcmp(number_format, orig_number_format), 'The file number storage formats are incosistant');
        %Update header
        if p.Results.append, headerSize = double([dims, szData(1:2), szData(3) + size(data,3)]);
        elseif p.Results.truncateToFrame, headerSize = double([dims, szData(1:2), p.Results.truncateToFrame]); %TODO this is not ideal in terms of the file is still gonna retain original size, we just never actually use the last parts
        else headerSize = double([dims,szData(1:3)]);
        end
        frewind(fid);
      else        
        headerSize = double([3, padarray(size(data),[0, 3-ndims(data)],1,'post')]);
      end
      fwrite(fid,headerSize,'double');
      
      %Write the header text (only on initial write %TODO, some asserts if it is correct)
      if ~to_modify
        headerStr = blanks(100);
        lPrefix = length(p.Results.prefix);
        headerStr(1:min(10,length(number_format))) = number_format(1:min(10,length(number_format)));
        headerStr(11:(10+min(lPrefix,70))) = p.Results.prefix(1:min(lPrefix,70));
        headerStr(end-length(p.Results.timestamp)+1:end) = p.Results.timestamp;

        fwrite(fid,headerStr,'char');
      end
      
      %Write the data
      if ~to_modify
        fwrite(fid,data(:),number_format);
      else %When we do modification to existing stack
        frewind(fid);
        %Go to first frame start
        dims = fread(fid, 1, 'double'); %get ndims
        szData = fread(fid,dims,'double')'; %get frame sizes
        fseek(fid,100,'cof'); %skip the header
        frameByteSkip = szData(1)*szData(2)*number_format_bytes;
        if p.Results.append %Continue writing from last frame
          %Skip all true frames
          fseek(fid,(szData(3)-size(data,3))*frameByteSkip,'cof');
          fwrite(fid,data(:),number_format);
        elseif p.Results.overwrite_frame
          %Skip enough frames
          fseek(fid,(p.Results.overwrite_frame-1)*frameByteSkip,'cof');
          %Write the given data from current position
          fwrite(fid,data(:),number_format);
        end
      end
      
      %Close the file
      fclose(fid);
      
end
