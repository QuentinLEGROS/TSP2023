function pprior=compute_pprior_LapMu(T,Nz,V,cl)


% Compute the log of the Lap prior
% 
% INPUT:
% T        : Current IF estimate
% Nx       : Number of admissible ridge position
% V        : indices of neighbors
% c        : Lap prior weight
%
% OUTPUT:
% pprior   : Log Lap prior
%
% Author: Q.Legros



[N,NN]=size(V);
Xm=200;
T1=ones(N,1)*(1:Nz); %N x Nz
pprior=zeros(N,Nz);
for t=1:NN
    X=((T1-T(V(:,t))*ones(1,Nz))).^2;
    X(isnan(X))=0;
    X(X>Xm)=Xm;
    pprior=pprior + X;
end
pprior=-2*cl*pprior;




