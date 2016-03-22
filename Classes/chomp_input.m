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


    function s = export_struct(obj)
      p = properties(obj);
      for i1 = 1:numel(p)
        s.(p{i1}) = obj.(p{i1});
      end
    end
    
  end
  
  methods (Static)
    %Check to set data pathes in a different environment on load
    function obj = loadobj(obj)
      %Look for enviroment variable
      addlistener(obj.opt, 'root_folder', 'PreSet',@(src, evnt)chomp_input.opt_change_pre(obj, src, evnt));
      addlistener(obj.opt, 'root_folder', 'PostSet',@(src, evnt)chomp_input.opt_change_post(obj, src, evnt));
      
      obj.opt.root_folder = getenv('CHOMP_ROOT_FOLDER'); %TODO make sure it is firing the event
    end
    
    %Change data pathes if options change
    function opt_change_pre(obj, src, evnt)
      if obj.opt.verbose > 1
        disp('opt_change_pre')
        disp(obj);
      end
      obj.data.raw_stack.Y.Source(1:(length(obj.opt.root_folder))) = []; %Just get the substring with the original root folder removed
      if isa(obj.data.proc_stack.Y, 'chomp_data')
        obj.data.proc_stack.Y.Source(1:(length(obj.opt.root_folder))) = [];
      end
    end
    
    function opt_change_post(obj, src, evnt)
      if obj.opt.verbose > 1
        disp('opt_change_post')
      end
      obj.data.raw_stack.Y.Source = [obj.opt.root_folder obj.data.raw_stack.Y.Source];
      if isa(obj.data.proc_stack.Y, 'chomp_data')
        obj.data.proc_stack.Y.Source = [obj.opt.root_folder obj.data.proc_stack.Y.Source];
      end
    end    
  end
  
end

