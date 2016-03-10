% SYMMYS - Last version of this code available at http://symmys.com/node/136
% Transforms raw moments into central moments (first central moment defined as expectation)

function mu=Raw2Central(mu_)



N=length(mu_);
mu=mu_;

for n=2:N
    mu(n) = ((-1)^n) * (mu_(1))^(n);
    for k=1:n-1
        mu(n) =  mu(n) + nchoosek(n,k) * ((-1)^(n-k)) * mu_(k) * (mu_(1))^(n-k);
    end
    mu(n) = mu(n)+ mu_(n);
end