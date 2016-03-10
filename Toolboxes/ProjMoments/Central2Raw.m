% SYMMYS - Last version of this code available at http://symmys.com/node/136
% Transforms central moments into raw moments (first central moment defined as expectation)

function mu_=Central2Raw(mu)

N=length(mu);
mu_=mu;

for n=2:N
    mu_(n) = ((-1)^(n+1)) * (mu(1))^(n);
    for k=1:n-1
        mu_(n) =  mu_(n) + nchoosek(n,k) * ((-1)^(n-k+1)) * mu_(k) * (mu_(1))^(n-k);
    end
    mu_(n) = mu_(n)+ mu(n);
end