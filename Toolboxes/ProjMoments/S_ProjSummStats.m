% SYMMYS - Last version of article and code available at http://symmys.com/node/136

% Project summary statistics to arbitrary horizons under i.i.d. assumption
% see Meucci, A. (2010) "Annualization and General Projection of Skewness, Kurtosis and All Summary Statistics"
% GARP Risk Professional, August, pp. 52-54 

clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=6;  % focus on first N standardized summary statistics
K=100; % projection horizon

% generate arbitrary distribution
J=100000;  % number of scenarios

a=-1; % shift
m=.2; % log-mean
s=.4; %log-sdev

Z=randn(J/2,1); 
Z=[Z;-Z]/std(Z,1);
X=a+exp(m+s*Z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute single-period standardized statistics and central moments
[ga,mu]=SummStats(X,N);

% compute single-period non-central moments
mu_=Central2Raw(mu);

% compute single-period cumulants
ka=Raw2Cumul(mu_);

% compute multi-period cumulants
Ka=K*ka;

% compute multi-period non-central moments
Mu_=Cumul2Raw(Ka);

% compute multi-period central moments
Mu=Raw2Central(Mu_);

% compute multi-period standardized statistics
Ga=Mu;
Ga(2)=sqrt(Mu(2));
for n=3:N
    Ga(n)=Mu(n)/(Ga(2)^n);
end