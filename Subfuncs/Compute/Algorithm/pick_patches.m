function [out] = pick_patches( datas, Hs, opts)
%PICK_PATCHES Returns patches to learn from, Y is the data tensor, H are the spatial locations, type is cell type number

%Make sure all the input are cell arrays (of possibly multiple datasets)
if ~iscell(datas), datas = {datas}; end
if ~iscell(Hs), Hs = {Hs}; end
if ~iscell(opts), opts = {opts}; end

py = cell(numel(datas),1);

parfor c1 = 1:numel(datas)
  data = datas{c1};
  H = Hs{c1};
  opt = opts{c1};
  py{c1} = cell(size(H,1),1);
  
  szY = chomp_size(data.proc_stack,'Y');

  patches = get_patch(data.proc_stack, opt, H);

  %Pick the patches
  for h1 = 1:size(H,1)
    curpy = get_n_order_patch(patches(:,:,:,h1), opt, szY);
    py{c1}{h1} = curpy;
  end


end

%Concatanate the results into a 1D array of cell arrays
out={};
for c1 = 1:numel(py)
out(end+1:end+numel(py{c1})) = py{1};

end
