function wc=Mod_Mstep_SEM(Y,Fc,wc,alpha,Ns,muc,P)


% Perform M-step
% 
% INPUT:
% Y           : Observation
% P           : Posterior distribution
% Fc          : data distribution
% wc          : current mixture weights
% alpha       : Dirichlet prior hyperparameter

% 
% OUTPUT:
% wc         : estimated mixture weights
%
% Author: Q.Legros

[N,~]=size(Y);
for n=1:N % Update timewise
    ind=find(P(n,:)>0.01*max(P(n,:))); % Select useful value of the posterior - speed-up the estimation
    p = P(n,ind)';
    
    wc(n,:) = Mod_NR_multi(Y(n,:)',wc(n,:),Fc,alpha,Ns,muc(n,:),p);
end