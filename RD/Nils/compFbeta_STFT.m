function [F_mat, Fbeta_mat, IntFbeta, Falpha_mat, varF]=compFbeta_STFT(beta,alpha,M,L,N)
%
% Compute the impulse response function (IRF) at the power beta. The Gabor
% kernel is approximated using a Gaussian distribution.
%
% INPUT:
% nwin             : length of the IRF
% meaF             : Gaussian mean
% sigmF            : Gaussian standard deviation
% beta             : hyperparameter
%
% OUTPUT:
% F_mat            : IRF (won't be use. Only needed for debuging)
% Fbeta_mat        : IRF at the power beta
% IntFbeta         : First term of the beta divergence

% Compute the sliding tfr window
val = transpose(Fh(-(M/4)+1:((3*M)/4), M, L ));

% Generate a 2D array for later use. Accelerate the computation of the
% cross entropy in 'online_2D'
F_mat = zeros(M,M/2);
F_mat(:,1) = val;
F_mat(:,1) = sqrt(F_mat(:,1));
for i = 1:(M/2)-1
    F_mat(:,i+1) = [F_mat(end,i);F_mat(1:end-1,i)];
end
F_mat = F_mat((M/4)+1:((3*M)/4),:); % Truncation to the same lenght than the data
F_mat = F_mat./sum(F_mat);% Normalization



MF=F_mat(:,M/4)'*(1:N)'; % mean for second moment computed below
varF=sqrt(F_mat(:,M/4)'*((1:N)'.^2)-MF.^2); % IRF Std for ridge removal

IntFbeta = sum(F_mat.^(1+beta)); % First term of the beta divergence
Fbeta_mat = F_mat.^beta; % Used in the second term of the beta divergence

Falpha_mat = (F_mat.^(1-alpha))+eps;% for Reyni divergence

