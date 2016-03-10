% SYMMYS - Last version of this code available at http://symmys.com/node/136
% Transforms cumulants into raw moments 

function mu_=Cumul2Raw(ka)

N=length(ka);
mu_=ka;

for n=1:N
    mu_(n) = ka(n);
    for k=1:n-1
        mu_(n) = mu_(n) + nchoosek(n-1,k-1) * ka(k)*mu_(n-k); 
    end
    
end