function [F_mat, Fbeta_mat, IntFbeta, Falpha_mat, varF, F_sAB]=compFbeta_STFT(beta,alpha,M,L,N)
%
% Compute the impulse response function (IRF) at the power beta. The Gabor
% kernel is approximated using a Gaussian distribution.
%
% INPUT:
% beta          : divergence hyperparameter
% alpha         : divergence hyperparameter
% M             : Number of frequential bin
% L             : analysis window size (in bin)
% N             : Signal length
% 
% 
% OUTPUT:
% F_mat         : Postulated obervation model
% Fbeta_mat     : Precomputed window for beta divergence
% IntFbeta      : Precomputed sum for beta divergence
% Falpha_mat    : Precomputed window for beta divergence
% varF          : Variance of the observation model

% Compute the sliding tfr window
val = transpose(Fh(-(M/4)+1:((3*M)/4), M, L ));

% Generate a 2D array for later use. Accelerate the computation of the
% cross entropy in 'online_2D'
F_mat = zeros(M,M/2);
F_mat(:,1) = val;
for i = 1:(M/2)-1
    F_mat(:,i+1) = [F_mat(end,i);F_mat(1:end-1,i)];
end
F_mat = F_mat((M/4):((3*M)/4)-1,:); % Truncation to the same lenght than the data
F_mat = F_mat + eps; % if beta or alpha < 0
F_mat = F_mat./sum(F_mat);% Normalization



MF=F_mat(:,M/4)'*(1:N)'; % mean for second moment computed below
varF=sqrt(F_mat(:,M/4)'*((1:N)'.^2)-MF.^2); % IRF Std for ridge removal

IntFbeta = sum(F_mat.^(1+beta)); % First term of the beta divergence
Fbeta_mat = F_mat.^beta; % Used in the second term of the beta divergence

Falpha_mat = (F_mat.^(1-alpha))+eps;% for Reyni divergence

F_sAB = (sum(F_mat.^(alpha+beta)));

