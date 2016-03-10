function out = raw2cum( inp, dim)
%RAW2CUM Summary of this function goes here
%   Detailed explanation goes here


%Put the moment dimension first
dims = 1:ndims(inp);
dims(1) = dim;
dims(dim) = 1;
inp = permute(inp, dims);

%Convert the raw moments into cumulants
szInp = size(inp);
inp = reshape(inp,szInp(1),[]);
out = inp;
for moms = 1:szInp(1)
  out(moms,:) = inp(moms,:);
  for i1=1:moms-1
      out(moms,:) = out(moms,:) - nchoosek(moms-1,i1-1) * out(i1,:)*inp(moms-i1,:); 
  end
end

%Permute the dimensions back
out = reshape(out,szInp);
out = ipermute(out,dims);
end

