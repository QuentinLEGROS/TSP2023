function Y = denSampling(X,P)

% Compute samples from P
% 
% INPUT:
% X        : Density support
% P        : Probability density function
%
% OUTPUT:
% Y        : samples
%
% Author: Q.Legros

[N,D]=size(P);
cum_P = cumsum(P,2);
z = rand(N,1);
A=(z*ones(1,D)<cum_P);
ind=D-sum(A,2)+1;
ind=min(max(ind,1),D);
Y=X(ind);
Y=Y(:);
