% SYMMYS - Last version of this code available at http://symmys.com/node/136
% Transforms raw moments into cumulants

function ka=Raw2Cumul(mu_)

N=length(mu_);
ka=mu_;

for n=1:N
    ka(n) = mu_(n);
    for k=1:n-1
        ka(n) = ka(n) - nchoosek(n-1,k-1) * ka(k)*mu_(n-k); 
    end
    
end