function [ ind1 ] = single_ind_shift( ind1, subs, inpsz, shifts )
%SINGLE_IND_SHIFT %ind 0 is a linear index in prod(inpsz)
% We want to compute ind1, a linear index in prod(inpsz) if we shift the
% dimensions according to shifts

%Check for every shift if we are still within the tensor, if not, just
%return NaN

%Check for all indices if they are still in the tensor after the shifts
subs = bsxfun(@plus, subs, shifts);
subs_toosmall = bsxfun(@ge, subs, ones(size(subs)));
subs_toolarge = bsxfun(@le, subs, repmat(inpsz,[size(subs,1),1]));
subs_good = min(subs_toosmall.*subs_toolarge,[],2);
subs_good(subs_good==0)=NaN;

ind1 = ind1.*subs_good;

ind1 = ind1 + shifts(1);
for dims = 2:length(shifts)
  ind1 = ind1 + shifts(dims)*prod(inpsz(1:(dims-1)));
end




end

