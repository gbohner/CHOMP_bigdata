function out = raw2cum_multivariate( inp )
%RAW2CUM Converts raw moments to cumulants for multivariate distributions for cell array inp, 
% where raw moments are stored within the individual cells inp{1} = mean,
% inp{2} = second moment etc

%Convert the raw moments into cumulants
out = inp;
for moms = 1:numel(inp)
  for i1=1:moms-1
      out{moms} = out{moms} - nchoosek(moms-1,i1-1) .* mply(out{i1}, inp{moms-i1}, 0); 
  end
end

end

