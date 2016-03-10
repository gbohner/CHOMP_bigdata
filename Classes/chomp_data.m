classdef chomp_data
  %CHOMP_DATA Store image stacks in a custom binary format with given
  %read-write functions
  
  properties
    Source %Source location of the binary file
  end
  
  methods %Constructor
    function obj = chomp_data(src, varargin)
      obj.Source = src;
      if nargin>1 %Construct with data input
        writeStack( src, varargin{1});
      end
    end
  end
  
  methods %Subsref and subassign
    function out = subsref(obj, subs)
      if strcmp(subs(1).type,'.')
        out = obj.(subs(1).subs);
      elseif strcmp(subs(1).type, '()'),'Wrong chomp_data subsref')
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
      end
    end
    
    function append(obj, data)
      writeStack(obj.Source,data,'append',1);
    end
    
    function szData = size(obj)
      fid = fopen(obj.Source,'r');
      dims = fread(fid, 1, 'double');
      szData = fread(fid,uint16(dims),'double')';
      fclose(fid);
    end
  end
  
end

