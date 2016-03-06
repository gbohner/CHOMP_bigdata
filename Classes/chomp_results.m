classdef chomp_results
  %CHOMP_OUTPUT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    opt %Chomp-options class object
    y %Mean processed image for visualization
    y_orig %Mean original image for visualization
    V %Variance processed image (we really don't need this though so explicitly?)
    UserMask %User defined mask to exclude areas we don't wanna do inference in
    H %Cell locations (index, not sub)
    X %Coefficients
    W %Basis functions %TODO: rename 
    L %Goodness of fit score
  end
  
  methods
    function obj = chomp_results(opt,W,H,X,L,y,y_orig,V)
      obj.opt = opt;
      obj.W = W;
      obj.H = H;
      obj.X = X;
      obj.L = L;
      obj.y = y;
      obj.y_orig = y_orig;
      obj.V = V;
    end
    
    function varargout = get_fields(obj,varargin)
      for i1 = 1:nargin-1
        varargout{i1} = obj.(varargin{i1});
      end
    end
  end
  
end

