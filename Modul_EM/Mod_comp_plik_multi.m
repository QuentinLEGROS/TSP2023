function C=Mod_comp_plik_multi(Y,Fc,wc,Ns)


% Compute the log likelihood 
% 
% INPUT:
% Y        : Observation
% Fc       : Convolution matrix of the Gaussian kernel - data distribution
% wc       : current mixture weights

% OUTPUT:
% C        : Log likelihood
%
% Author: Q.Legros


[N,M]=size(Y);
[~,Nx,~]=size(Fc);
C=zeros(N,Nx);
for n=1:N
    A =(1-sum(wc(n,:)))/M;
    for ns = 1:Ns
        A = A+(wc(n,ns)*Fc);
    end
    C(n,:)=sum((Y(n,:)'*ones(1,Nx)).*log(A+eps),1);

end




