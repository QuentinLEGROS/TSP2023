function wc=Mod_Mstep_multi(Y,Fc,wc,alpha,Ns,muc)


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

[M,Nx]=size(Fc);
[N,~]=size(Y);
for n=1:N % Update pixelwise
    if sum(Y(n,:))==0 % If spectrogram is zeros value : max the posterior relies to max the prior
        wc(n,:)=(alpha(1:Ns)-1)/(sum(alpha)-length(alpha)); % mode of the Dirichlet distribution
    else % else : max the posterior relies to max the prior
        
    wc(n,:) = Mod_NR_multi(Y(n,:)',wc(n,:),Fc,alpha,Ns,muc(n,:));
    end
end