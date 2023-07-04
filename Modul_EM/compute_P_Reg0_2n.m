function [P_out,Mu_out]=compute_P_Reg0_2n(plik0,Mu,c,reg_mu,N_samp,Ns,step_r,step_v)

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
ind0=cell(2,1);

for ns = 1:Ns
    for t=1:N_samp
        ind0=cell(2,1);
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
                ppost=plik+pprior;
                ppost=ppost-max(ppost,[],2)*ones(1,Nx);
                post=exp(ppost);
                ppost=post./(sum(post,2)*ones(1,Nx));
                Mu2(ind0{r},ns)=denSampling(1:Nx,ppost);
        end
    end
%     Mu2 = mean(Mu2(3:end-2,:,:),3);
    
    
    for r=randperm(2)
    ind0{r}=find(Z==r-1);
    if strcmp(reg_mu,'TV')
        V=[ind0{r}-1,ind0{r}+1];
        pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
    elseif strcmp(reg_mu,'Lap')
        V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
        pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
    end
    plik=P1(ind0{r},:);
     pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
    C((ind0{r}),:)=plik+pprior;
    end
end



C2=C(3:end-2,:);
C2=reshape(C2,N,Nx);
C2=C2-max(C2,[],2)*ones(1,Nx);
P=exp(C2);
P=P./(sum(P,2)*ones(1,Nx));
P_out = P;



% for t=1:N_samp
%     for ns = 1:Ns
%         ind0=cell(2,1);
%         C = zeros(N+4,Nx);
%         for r=randperm(2)
%             ind0{r}=find(Z==r-1);
%             plik=P1(ind0{r},:);
%             if strcmp(reg_mu,'TV')
%                 V=[ind0{r}-1,ind0{r}+1];
%                 pprior= compute_pprior_TVMu(Mu2(:,ns),Nx,V,c);
%             elseif strcmp(reg_mu,'Lap')
%                 V=[ind0{r}-1,ind0{r}-2,ind0{r}+1,ind0{r}+2];
%                 pprior= compute_pprior_LapMu(Mu2(:,ns),Nx,V,c);
%             end
%         C((ind0{r}),:)=plik+pprior;
%         cprior((ind0{r}),:) = pprior;
%         cplik((ind0{r}),:) = plik;
%     end
% 
%     C2=C(3:end-2,:);
%     C2=reshape(C2,N,Nx);
%     C2=C2-max(C2,[],2)*ones(1,Nx);
%     P=exp(C2);
%     P=P./(sum(P,2)*ones(1,Nx));
%     P_out = P_out+P;
% 
%     end
%     P_out=P_out./(sum(P_out,2)*ones(1,Nx));
%     Mu2 = compute_itSMAP(P_out,Ns,step_r,step_v);
%     Mu2 = [NaN.*ones(2,2);Mu2;NaN.*ones(2,2)];
% end
% 
% %     C = ppost;
%     C2=C(3:end-2,:);
%     C2=reshape(C2,N,Nx);
%     C2=C2-max(C2,[],2)*ones(1,Nx);
%     P=exp(C2);
%     P=P./(sum(P,2)*ones(1,Nx));
%     P_out = P_out+P;

    
P_out=P_out./(sum(P_out,2)*ones(1,Nx));
Mu_out = compute_itSMAP(P_out,Ns,step_r,step_v);

Mu_out = Mu2(3:end-2,:);

% figure
% subplot(2,1,1)
% plot(Mu2)
% subplot(2,1,2)
% plot(Mu_out)

[~,nn]=sort(mean(Mu_out,1));
Mu_out = Mu_out(:,flip(nn));

% disp('la')

