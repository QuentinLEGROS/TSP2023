function [mcav_ext,s2cav_ext]=online_2D(Y,LFM,Fbeta,Falpha,ds,mcav_ext,s2cav_ext,IntFbeta,beta,alpha,div)
%
% Online approach. Pseudo-Bayesian method for estimating the current (for 
% the current time considered) position of the ridge. The pseudo-likelihood
% term is the beta divergence and the temporal prior a Gaussian pdf
% centerer around the previous mean estimated (at time t-1) and with a
% time-increasing standard deviation (Gaussian random walk prior).
%
% INPUT:
% Y              : data
% tt             : support
% Fbeta          : IRF at power beta
% ds             : variance of the Gaussian random walk
% mcav_ext       : mean of the prior
% s2cav_ext      : standard deviation of the prior
% 
% OUTPUT:
% mcav_ext       : estimated mean (position of the ridge by MMSE)
% s2cav_ext      : estimated variance of the estimate

T = length(Y);
%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT :  need to add 'eps' value to F_mat to compute the log since it
% contains zeros values.

%% Compute log-divergence (pseudo-likelihood)
if div == 1 % KL divergence
    U = (Y*LFM); % Cross KL term
%     U = U - sum(Y.*log(Y)); % Add expectation on log data
elseif div == 2 % beta divergence
    U=(Y*Fbeta); % beta cross entropy
    U = ((1+beta)/beta)*U - IntFbeta;% - ((1/beta)*Y.^(1+beta));
elseif div == 3 % Renyi divergence
    U = -log((Y.^alpha)*(Falpha))./(alpha-1);
else % KL divergence
    U = (Y*LFM); % Cross KL term
%     U = U - sum(Y.*log(Y)); % Add expectation on log data
end


% figure(40)
% subplot(1,2,1)
% hold on
% plot(Y)
% plot(U)
% plot(U2)
% plot(U3)
% hold off
% legend('data','KL','beta','alpha')
% 




%% Compute prior 
SS=s2cav_ext+ds;
F1=normpdf(1:T,mcav_ext,sqrt(SS));
%% Compute posterior (combine pseudo-likelihood and prior)
U=U+log(F1+eps);
U=U-max(U,[],2)*ones(1,T);
P=exp(U);
P=P./sum(P);

% U2=U2+log(F1+eps);
% U2=U2-max(U2,[],2)*ones(1,T);
% P2=exp(U2);
% P2=P2./sum(P2);
% 
% U3=U3+log(F1+eps);
% U3=U3-max(U3,[],2)*ones(1,T);
% P3=exp(U3);
% P3=P3./sum(P3);
% 
% 
% subplot(1,2,2)
% hold on
% plot(Y)
% plot(P)
% plot(P2)
% plot(P3)
% hold off
% legend('data','KL','beta','alpha')

%% Compute the mean and variance of the posterior
% MAP
% [~,mcav_ext] = max(P);

% MMSE
mcav_ext=P*(1:T)';
s2cav_ext=P*((1:T)'.^2)-mcav_ext.^2;


