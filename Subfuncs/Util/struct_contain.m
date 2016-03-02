function output = struct_contain(s1, s2, varargin)
%STRUCT_CONTAIN Returns true if s2 has all fields in s1 with equal values.
%   as a 3rd one can set a cell array of specific field names to check for
%   Furthermore, as a 4th argument, you can specify specific field names to ignore

  %Return true by default, we check if anything violates it
  output = 1;

  %Select which fields to check for
  if nargin==2
    names = fieldnames(s1);
  elseif nargin==3
    names = varargin{1};
  elseif nargin==4
    names = fieldnames(s1);
    for i1=1:length(varargin{2})
      ind=find(ismember(names,varargin{2}{i1}));
      if ~isempty(ind)
        names(ind) = [];
      end
    end
  end

  %Check if s2 indeed has all field names
  names2 = fieldnames(s2);
  for i1=1:length(names)
    ind=find(ismember(names2,names{i1}));
    if isempty(ind)
      %Return 0 if a field_name in s1 is not present in s2
      output=0;
      return;
    end
  end

  %Check if all the values are equal
  for i1=1:length(names)
    if ~isequal(s1.(names{i1}), s2.(names{i1}))
      %Check if isequal doesn't work because it is function handles
      if isa(s1.(names{i1}),'function_handle') && isa(s2.(names{i1}),'function_handle')
        if isequal(func2str(s1.(names{i1})), func2str(s2.(names{i1}))), continue; end
      end
      output=0;
      return
    end
  end

end

