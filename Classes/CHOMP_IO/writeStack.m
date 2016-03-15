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
      
      number_format = 'uint16';
      tmp = 1; tmp = cast(tmp,number_format); tmp = whos('tmp');
      number_format_bytes = tmp.bytes;
      
      %File format:
      p = inputParser();
      p.addRequired('path',@ischar)
      p.addRequired('data',@isnumeric)
      p.addParameter('prefix','DEFAULT_PREFIX',@ischar)
      p.addParameter('timestamp',datestr(now, 30), @ischar)
      p.addParameter('append',0,@(x)all([isnumeric(x),exist(path,'file')]));
      p.addParameter('overwrite_frame',0,@(x)all([isnumeric(x),exist(path,'file')])); %overwrite frame #overwrite_frame
      p.parse(path, data,varargin{:})
      
      %Open the file at the given location
      if p.Results.append || p.Results.overwrite_frame
        fid = fopen(p.Results.path,'rb+');
      else
        fid = fopen(p.Results.path,'w');
      end
      
      %Convert the data to required number format
      p.Results.data = cast(p.Results.data,number_format);
      
      %Write into the header [ndims(data), size(data)]
      if p.Results.append
        %Get the size of already existing data
        frewind(fid);
        dims = fread(fid, 1, 'double');
        szData = fread(fid,uint16(dims),'double')';
        assert(all(szData(1:2)==[size(p.Results.data,1),size(p.Results.data,2)]),'Dimension for appending are inconsistant');
        headerSize = double([dims, szData(1:2), szData(3) + size(p.Results.data,3)]);
        frewind(fid);
      else        
        headerSize = double([3, padarray(size(p.Results.data),[0, 3-ndims(p.Results.data)],1,'post')]);
      end
      fwrite(fid,headerSize,'double');
      
      %Write the header text (unless append)
      if ~p.Results.append
        headerStr = blanks(100);
        lPrefix = length(p.Results.prefix);
        headerStr(1:min(10,length(number_format))) = number_format(1:min(10,length(number_format)));
        headerStr(11:(10+min(lPrefix,70))) = p.Results.prefix(1:min(lPrefix,70));
        headerStr(end-length(p.Results.timestamp)+1:end) = p.Results.timestamp;

        fwrite(fid,headerStr,'char');
      end
      
      %Write the data
      if ~p.Results.append
        fwrite(fid,p.Results.data(:),number_format);
        
      elseif p.Results.append %Continue writing from eof
        fseek(fid,0,'eof');
        fwrite(fid,p.Results.data(:),number_format);
      elseif p.Results.overwrite_frame
        frewind(fid);
        %Go to first frame start
        dims = fread(fid, 1, number_format); %get ndims
        szData = fread(fid,uint16(dims),number_format)'; %get frame sizes
        fseek(fid,100,'cof'); %skip the header
        %Skip enough frames
        frameByteSkip = szData(1)*szData(2)*number_format_bytes;
        fseek(fid,uint16((p.results.overwrite_frame-1)*frameByteSkip),'cof');
        %Write the given data from current position
        fwrite(fid,p.Results.data(:),number_format);
      end
      
      %Close the file
      fclose(fid);
      
end
