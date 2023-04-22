function [P_out,Mu_out]=compute_P_Reg2_2n(plik0,Mu,c,reg_mu,N_samp,Ns,step_r,step_v)

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
%     for t=1:N_samp
%         for r=randperm(2)
%             ind0{r}=find(Z==r-1); % chosen even or odd elements
%             plik=P1(ind0{r},:); % likelihood values in the chosen bins
%             if strcmp(reg_mu,'TV') % if TV
%                 V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
%                 pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
%             elseif strcmp(reg_mu,'Lap') % if Lap
%                 V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
% %                 V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
%                 pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
%             end
%             ppost=plik+pprior;
%             ppost=ppost-max(ppost,[],2)*ones(1,Nx);
%             post=exp(ppost);
%             ppost=post./(sum(post,2)*ones(1,Nx));
%             Mu2(ind0{r},ns)=denSampling(1:Nx,ppost);
%         end
%     end
    
    

    C=zeros(N+4,Nx);
    %% Compute posterior from samples after convergence of the priors
    for r=randperm(2)
        ind0{r}=find(Z==r-1);
        plik=P1(ind0{r},:);
            if strcmp(reg_mu,'TV')
                V=[ind0{r}-1,ind0{r}+1];
                pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
            elseif strcmp(reg_mu,'Lap')
                V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
%               V=[ind0{r}-1,ind0{r}+1]; % neighborhoods of the chosen bins
                pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
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
    P_out = P_out+P;
    
%     figure
%     imagesc(P)
    
% [~,Mu_out(:,ns)] = max(P,[],2);


end
% figure(1)
% imagesc(P_out)
% pause(0.01)
P_out=P_out./(sum(P_out,2)*ones(1,Nx));

% Mu_out = compute_SMAP(P_out,Ns,Mu);
Mu_out = compute_itSMAP(P_out,Ns,step_r,step_v);


[~,nn]=sort(mean(Mu_out,1));
Mu_out = Mu_out(:,flip(nn));

% disp('la')




