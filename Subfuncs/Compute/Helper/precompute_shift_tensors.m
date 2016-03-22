function GPT = precompute_shift_tensors(opt)

savefile = get_path(opt, 'precomputed');

% Compute shift tensors for all moments and shifts
GPT = cell(2*opt.m-1, 2*opt.m-1, opt.mom);

if opt.verbose
  fprintf('\nComputing the interaction tensor...\n');
end

for mom1 = 1:opt.mom
  inpuf = [ones(1,2*mom1)*opt.m]; %unfolded dimensions in non-feature-space
  outDim = opt.m^(2*mom1);
  
  charcount = 0;
  %Iterate through shifts
  for s1 = 1:(2*opt.m-1)
    if opt.verbose > 1
      if charcount>0, for c1 = 1:charcount, fprintf('\b'); end; end
      charcount = fprintf('Progress is currently %d/%d moment, %d/%d shift\n',mom1, opt.mom, s1, 2*opt.m-1);
    end
 
    curm = opt.m; %just to use in the parfor not having to pass through the whole opt struct
    for s2 = 1:(2*opt.m-1)
      %unfold, shift corresponding dimensions, fold again, then multiply
      %Compute the vectors of shifting
        s1_vec = zeros(1,length(inpuf));
        s2_vec = zeros(1,length(inpuf));
        s1_vec(1:2:(2+2*(mom1-1))) = s1-curm;
        s2_vec(2:2:(2+2*(mom1-1))) = s2-curm;
        shift_vec = s1_vec + s2_vec;
        oldInd = (1:outDim)';
        [tmp{1:length(inpuf)}] = ind2sub(inpuf, oldInd); %Get the sub-description for all linear indices
        subs = cell2mat(tmp);
        newInd = single_ind_shift(oldInd, subs,inpuf, shift_vec);

        oldInd(isnan(newInd)) = [];
        newInd(isnan(newInd)) = [];
        
        
      GPT{s1,s2,mom1} = [oldInd, newInd];
      % Shift every dimension with s1-s2-s1-s2-etc then multiply with the
      % inverse
    end
  end
end

save(savefile, 'GPT', '-v7.3');