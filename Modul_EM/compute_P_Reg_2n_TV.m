function [P_out,Mu_out]=compute_P_Reg_2n_TV(plik0,Mu,c,Ns,step_r,step_v)

% Estimate the posterior distribution using samples
% 
% INPUT:
% plik0     : Log likelihood
% T         : Current depth estimate
% c         : prior weight
%
% OUTPUT:
% P         : Posterior distribution
%
% Author: Q.Legros


%% Init
[N,Nx]=size(plik0);
Z=NaN*ones(N+4,1);
P1=zeros(N+4,Nx);
P1(3:end-2,:)=plik0;
Mu1=NaN*ones(N+4,Ns);
Mu1(3:end-2,:)=Mu;
Z(3:2:end-2)=0;
Z(4:2:end-2)=1;
P_out = zeros(size(plik0));
pprior = zeros(Nx,Nx);
%% Compute posterior and sample
Mu2=Mu1;
for ns = 1:Ns
    ind0=cell(2,1);
    C=zeros(N+4,Nx);
    %% Compute posterior from samples after convergence of the priors
    for r=randperm(2)
        ind0{r}=find(Z==r-1);
        plik=P1(ind0{r},:);
        V=[ind0{r}-1,ind0{r}+1];
        pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
        C((ind0{r}),:)=plik+pprior;
    end

    C=C(3:end-2,:);
    C=reshape(C,N,Nx);
    C=C-max(C,[],2)*ones(1,Nx);
    P=exp(C);
    P=P./(sum(P,2)*ones(1,Nx));
    P_out = P_out+P;
    
end
P_out=P_out./(sum(P_out,2)*ones(1,Nx));

Mu_out = compute_itSMAP(P_out,Ns,step_r,step_v);

[~,nn]=sort(mean(Mu_out,1));
Mu_out = Mu_out(:,flip(nn));




