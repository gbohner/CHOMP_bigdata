classdef chomp_input < handle
  %CHOMP_INPUT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetObservable)
    opt %Class that stores all the options
  end
  
  properties
    data %Struct that stores all the data or pathes to the virtual stacks containing the data
    y %Mean processed image for visualization
    y_orig %Mean original image for visualization
    V %Variance processed image (we really don't need this though so explicitly?)
    UserMask %User defined mask to exclude areas we don't wanna do inference in
  end
  
  methods
    function obj = chomp_input(opt1,data1,y1,y_orig1,V1)
      obj.opt = opt1;
      obj.data = data1;
      obj.y = y1;
      obj.y_orig = y_orig1;
      obj.V = V1;
      
      %Change the pathes when options change
      addlistener(obj.opt, 'root_folder', 'PreSet',@(src, evnt)chomp_input.opt_change_pre(obj, src, evnt));
      addlistener(obj.opt, 'root_folder', 'PostSet',@(src, evnt)chomp_input.opt_change_post(obj, src, evnt));
    end

    %Check to set data pathes in a different environment on load
    function loadobj()
      %Look for enviroment variable
      obj.opt.root_folder = getenv('CHOMP_ROOT_FOLDER'); %TODO make sure it is firing the event
    end
    
    function s = export_struct(obj)
      p = properties(obj);
      for i1 = 1:numel(p)
        s.(p{i1}) = obj.(p{i1});
      end
    end
    
  end
  
  methods (Static)
    %Change data pathes if options change
    function opt_change_pre(obj, src, evnt)
      disp('opt_change_pre')
      disp(obj);
      obj.data.raw_stack = obj.data.raw_stack.Properties.Source((length(obj.opt.root_folder)+1):end); %Just get the substring with the original root folder removed
      if isa(obj.data.proc_stack, 'matlab.io.MatFile')
        obj.data.proc_stack = obj.data.proc_stack.Properties.Source((length(obj.opt.root_folder)+1):end); %Just get the substring with the original root folder removed
      end
    end
    
    function opt_change_post(obj, src, evnt)
      disp('opt_change_post')
      obj.data.raw_stack = matfile([obj.opt.root_folder obj.data.raw_stack]);
      if isa(obj.data.proc_stack, 'char')
        obj.data.proc_stack.Properties.Source = matfile([obj.opt.root_folder obj.data.proc_stack]);
      end
    end    
  end
  
end

