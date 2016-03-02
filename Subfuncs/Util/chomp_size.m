function varargout = chomp_size( obj, varargin )
%SIZE Summary of this function goes here
%   Detailed explanation goes here

if nargin >1 
  field = varargin{1};
end

if isa(obj,'matlab.io.MatFile')
  varargout{:} = size(obj, field);
elseif isa(obj,'struct')
  varargout{:} = size(obj.(field));
elseif isa(obj,'matlab.io.datastore.ImageDatastore')
  if numel(obj.Files) == 1
    varargout{:} = size(readimage(obj,1)); %Single tif stack
  else
    varargout{:} = size(readimage(obj,1));
    varargout{1}(end+1) = numel(obj.Files); %How many frames stored
  end
end

end

