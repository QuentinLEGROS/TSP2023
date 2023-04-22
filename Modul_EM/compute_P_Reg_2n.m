function [P,Mu_out]=compute_P_Reg_2n(plik0,Mu,c,reg_mu,N_samp,Ns)

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
C=zeros(N+4,Nx);
pprior = zeros(Nx,Nx);
%% Compute posterior and sample
Mu2=Mu1;
for ns = 1:Ns
    ind0=cell(2,1);
    for t=1:N_samp
        for r=randperm(2)
            ind0{r}=find(Z==r-1); % chosen even or odd elements
            plik=P1(ind0{r},:); % likelihood values in the chosen bins
            if strcmp(reg_mu,'TV') % if TV
                V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
                pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
            elseif strcmp(reg_mu,'Lap') % if Lap
%                 V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
                V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
                pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
            end
            ppost=plik+pprior;
            ppost=ppost-max(ppost,[],2)*ones(1,Nx);
            post=exp(ppost);
            ppost=post./(sum(post,2)*ones(1,Nx));
            Mu2(ind0{r},ns)=denSampling(1:Nx,ppost);
        end
    end
end
% figure;subplot(3,1,1);plot(plik(40,:));subplot(3,1,2);plot(pprior(40,:));subplot(3,1,3);plot(ppost(40,:));

pprior = zeros(Nx,Nx);
%% Compute posterior from samples after convergence of the priors
% Mu_out=Mu2;
for r=randperm(2)
    ind0{r}=find(Z==r-1);
    plik=P1(ind0{r},:);
        if strcmp(reg_mu,'TV')
            V=[ind0{r}-1,ind0{r}+1];
            for ns = 1:Ns
                pprior=pprior + compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
            end        
        elseif strcmp(reg_mu,'Lap')
              V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
%             V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
            for ns = 1:Ns
                pprior=pprior + compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
            end        
        end
    C((ind0{r}),:)=plik+pprior;
    cprior((ind0{r}),:) = pprior;
    cplik((ind0{r}),:) = plik;
end

C=C(3:end-2,:);
C=reshape(C,N,Nx);
C=C-max(C,[],2)*ones(1,Nx);
P=exp(C);
P=P./(sum(P,2)*ones(1,Nx));
% for nn=1:N
% figure(1);subplot(3,1,1);plot(cplik(nn,:));subplot(3,1,2);plot(cprior(nn,:));subplot(3,1,3);plot(P(nn,:));
% pause(0.1)
% end
% disp('ok')

Mu_out=Mu2(3:end-2,:);
Mu_out2 = compute_SMAP(P,Ns);

figure(5)
subplot(2,2,1)
plot(Mu_out(:,1))
subplot(2,2,2)
plot(Mu_out(:,2))
subplot(2,2,3)
plot(Mu_out2(:,1))
subplot(2,2,4)
plot(Mu_out2(:,2))
disp('la')
% Mu_out=Mu_out(:);





