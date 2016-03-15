classdef chomp_data
  %CHOMP_DATA Store image stacks in a custom binary format with given
  %read-write functions
  
  properties
    Source %Source location of the binary file
  end
  
  methods %Constructor
    function obj = chomp_data(src, varargin)
      obj.Source = src;
      if nargin>1 %Pass the rest of the arguments on to writeStack
          writeStack( src, varargin{:});
      end
      
    end
  end
  
  methods %end, Subsref and subassign
    function out = end(obj,k,~)
      %Returns the size of obj for the kth dimension
      szData = size(obj);
      out = szData(k);
    end
    
    function out = subsref(obj, subs)
      if strcmp(subs(1).type,'.')
        out = obj.(subs(1).subs);
      elseif strcmp(subs(1).type, '()')
        assert(numel(subs(1).subs)==numel(size(obj)),'Wrong chomp_data subsref #dims');
        szData = size(obj);
        if strcmp(subs(1).subs{3},':')
          frames = 1:szData(3);
        else
          frames = subs(1).subs{3};
        end
        if strcmp(subs(1).subs{1},':') && strcmp(subs(1).subs{2},':')
          %Getting full frames
          out = readStack(obj.Source,frames);
        else
          if strcmp(subs(1).subs{1},':')
            patch.x = 1:szData(1);
          else
            patch.x = subs(1).subs{1};
          end
          if strcmp(subs(1).subs{2},':')
            patch.y = 1:szData(2);
          else
            patch.y = subs(1).subs{2};
          end
          out = readStack(obj.Source,frames,'patch',patch);
        end
      else
        error(sprintf('The type of subsref you tried to use, %s, is not handled by chomp_data',subs(1).type));
      end
    end
    
    function obj = subsasgn(obj,subs,val)
      if strcmp(subs(1).type,'.')
        obj.(subs(1).subs) = val;
      elseif strcmp(subs(1).type,'()')
        %Only handle overwriting a set of frames, appending or truncating
        assert(strcmp(subs(1).subs{1},':') && strcmp(subs(1).subs{2},':'), 'Only full frames can be overwritten, patches cannot')
        szData = size(obj);
        new_frames = subs(1).subs{3};
        if new_frames(1) == szData(3)+1
          %Check for appending
          assert(all(new_frames == new_frames(1):new_frames(end)),'Cannot append sparse frames');
          writeStack(obj.Source,val,'append',1,'number_format',obj.numFormat);
        elseif all(new_frames <= szData(3)) && ~isempty(val)
          %Check for overwriting
          assert(numel(new_frames) == size(val,3),'Attempting to overwrite with incorrect amount of data')
          for t1 = 1:numel(new_frames)
            writeStack(obj.Source,val(:,:,t1),'overwrite_frame',new_frames(t1),'number_format',obj.numFormat);
          end
        elseif (new_frames(end) == szData(3)) && isempty(val)
          %Check for truncating
          assert(all(new_frames == new_frames(1):new_frames(end)),'Cannot truncate sparse frames');
          assert(new_frames(1)>1, 'Cannot truncate all of the data');
          writeStack(obj.Source,[],'truncateToFrame',new_frames(1)-1);
        else
          error('The type of subsasgn you tried to use is not handled by chomp_data');
        end
      else
        error(sprintf('The type of subsasgn you tried to use, %s, is not handled by chomp_data',subs(1).type));
      end      
    end
    
    function overwrite_frame(obj, data, frameNum, varargin)
      writeStack(obj.Source,data,'overwrite_frame',frameNum, varargin{:});
    end
    
    function append(obj, data, varargin)
      writeStack(obj.Source,data,'append',1, varargin{:});
    end
    
    function out = size(obj)
      fid = fopen(obj.Source,'r');
      dims = fread(fid, 1, 'double');
      out = fread(fid,uint16(dims),'double')';
      fclose(fid);
    end
    
    function out = numFormat(obj)
      fid = fopen(obj.Source,'r');
      dims = fread(fid, 1, 'double');
      fseek(fid,uint16(dims*8),'cof'); %skip dimensions
      out = strtrim(char(fread(fid,10,'char')'));
      fclose(fid);
    end
  end
  
end

