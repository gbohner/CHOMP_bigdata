function writeStack( path, data, varargin )
      %Writes the data into binary file(s)
      %File format:
      p = inputParser();
      p.addRequired('path',@ischar)
      p.addRequired('data',@isfloat)
      p.addParameter('prefix','TEST_PREFIX_TOO_LONG',@ischar)
      p.addParameter('timestamp',datestr(now, 30), @ischar)
      p.addParameter('append',0,@(x)all([isnumeric(x),exist(path,'file')]));
      p.parse(path, data,varargin{:})
      
      %Open the file at the given location
      if p.Results.append
        fid = fopen(p.Results.path,'rb+');
      else
        fid = fopen(p.Results.path,'w');
      end
      
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
        headerStr(1:min(lPrefix,70)) = p.Results.prefix(1:min(lPrefix,70));
        headerStr(end-length(p.Results.timestamp)+1:end) = p.Results.timestamp;

        fwrite(fid,headerStr,'char');
      end
      
      %Write the data
      if ~p.Results.append
        fwrite(fid,p.Results.data(:),'double');
        
      else %Continue writing from eof
        fseek(fid,0,'eof');
        fwrite(fid,p.Results.data(:),'double');
      end
      %Close the file
      fclose(fid);
      
end
