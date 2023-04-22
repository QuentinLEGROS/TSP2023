function [P,Mu_out]=compute_P_Lap_2n(plik0,Mu,c,reg_mu,N_samp)

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
Mu1=NaN*ones(N+4,1);
Mu1(3:end-2,1)=Mu;
Z(3:2:end-2)=0;
Z(4:2:end-2)=1;
C=zeros(N+4,Nx);

%% Compute posterior and sample
Mu2=Mu1;
ind0=cell(2,1);
for t=1:N_samp
    for r=randperm(2)
        ind0{r}=find(Z==r-1);
        plik=P1(ind0{r},:);
        if strcmp(reg_mu,'TV')
            V=[ind0{r}-1,ind0{r}+1];
            pprior=compute_pprior_TVMu(Mu2(:),Nx,V,c); % la
        elseif strcmp(reg_mu,'Lap')
%             V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
            V=[ind0{r}-1,ind0{r}+1];
            pprior=compute_pprior_LapMu(Mu2(:),Nx,V,c); % la
        end
        ppost=plik+pprior;
        ppost=ppost-max(ppost,[],2)*ones(1,Nx);
        post=exp(ppost);
        ppost=post./(sum(post,2)*ones(1,Nx));
        Mu2(ind0{r})=denSampling(1:Nx,ppost);
    end
end
% C = ppost(2:end-1,:);
pprior = zeros(size(pprior));
%% Compute posterior from samples
Mu_out=Mu2;
for r=randperm(2)
    ind0{r}=find(Z==r-1);
    plik=P1(ind0{r},:);
        if strcmp(reg_mu,'TV')
            V=[ind0{r}-1,ind0{r}+1];
            pprior=compute_pprior_TVMu(Mu2(:),Nx,V,c); % la
        elseif strcmp(reg_mu,'Lap')
            V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
            pprior=compute_pprior_LapMu(Mu2(:),Nx,V,c); % la
        end
    C((ind0{r}),:)=plik+pprior;
end

C=C(3:end-2,:);
C=reshape(C,N,Nx);
C=C-max(C,[],2)*ones(1,Nx);
P=exp(C);
P=P./(sum(P,2)*ones(1,Nx));
Mu_out=Mu_out(3:end-2,1);
Mu_out=Mu_out(:);


