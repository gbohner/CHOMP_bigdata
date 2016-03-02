function output = get_path( opt,varargin)
%GET_PATH Gets the path for certain important files, by default the intermediate input data file, with different secondary options you can modify that

if nargin==1
  output = [opt.input_folder opt.file_prefix '_' opt.timestamp '.mat'];
else
  switch varargin{1}
    case 'output_iter'
      %Require which output iteration as varargin{2}
      output = [opt.output_folder opt.file_prefix '_' opt.timestamp '_iter_' num2str(varargin{2}) '.mat'];
    case 'precomputed'
      output = [opt.precomputed_folder 'Precomputed_shift_tensor_window_' num2str(opt.m) '_moment_' num2str(opt.mom) '.mat'];
    case 'virtual_stack'
      output = [opt.input_folder opt.file_prefix '_' opt.timestamp '_virtual_stack.mat'];
    case 'raw_virtual_stack'
      output = [opt.input_folder opt.file_prefix '_virtual_stack_raw.mat'];
    otherwise
       error('CHOMP:util:nopath', 'The type of path you want get_path() to return is not implemented');
  end
end
      


end

