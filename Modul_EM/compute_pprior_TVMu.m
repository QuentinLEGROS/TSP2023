function pprior=compute_pprior_TVMu(T,Nx,V,c)


% Compute the log of the TV prior
% 
% INPUT:
% T        : Current IF estimate
% Nx       : Number of admissible ridge position
% V        : indices of neighbors
% c        : TV prior weight
%
% OUTPUT:
% pprior   : Log TV prior
%
% Author: Q.Legros



[N,NN]=size(V);
Xm=200;
T1=ones(N,1)*(1:Nx); %N x Nz
pprior=zeros(N,Nx);
for t=1:NN
    X=abs(T1-T(V(:,t))*ones(1,Nx));
    X(isnan(X))=0;
    X(X>Xm)=Xm;
    pprior=pprior + X;
end

pprior=-2*c*pprior;



