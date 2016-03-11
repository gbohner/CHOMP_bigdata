function out = lmsvd_funcs( x, params )
%LMSVD_FUNCS Multiples x directly with data extracted (possibly from multiple locations),
% only storing the results in the memory

%TODO: LMSVD joint learning, but: inconsistent number of samples / dataset

%Cell arrays containing all the required data
[stacks, opts, Hs, func_to_use] = params{:};

results = cell(numel(stacks),1);

parfor n1 = 1:numel(stacks)
  patches = get_patch(stacks{n1}, opts{n1}, Hs{n1});
  szP = size(patches);
  results{n1} = feval(funcs, reshape(patches,prod(szP(1:2)),[]), x);
end

  function z = lmply(A,x)
    z = A*x;
  end

  function z = trans_mply(A,x)
    z = A'*x;
  end



end

