function [A_hat] = Oracle_Amp_DF(tfr,Ncomp,M,L,tf0)
%
% Main algorithm: estimate the ridge position and the variance of the
% posterior distribution to propagate information to the next time sequence
%
% INPUT:
% tfr           : Time-frequency represerntation of the MCS
% Ncomp         : Number of component
% M             : Number of frequential bin
% L             : analysis window size (in bin)
% tf0           : Oracle instantaneous frequency
%
%
% OUTPUT:
% A_hat      : Estimated amplitude


if ~exist('ds', 'var')
 ds=3;% variance of the random walk in the temporal model
end
if ~exist('Pnei', 'var')
 Pnei=10;% default KL
end


[M,N] = size(tfr); % Extract dimenssions

x_hat = zeros(Ncomp,N);
A_hat = zeros(N,Ncomp);

%% compute mask
Mask_out = compMask(tf0,Pnei,M/2,0);

% med_win = [-2 -1 0 1 2];
% A_hat_temp = zeros(1,length(med_win));

%% Estimation amplitude - DF code
for Nc = 1:Ncomp
    x_hat(Nc,:) = real(rectfrgab(tfr .* Mask_out(:,:,Nc), L, M));
    for n = 1:N 
%         for indw = 1:length(med_win)
%             tf1 = min(max(tf0(n,Nc)+med_win(indw),1),N);
            A_hat(n,Nc) = Compute_A_DF2(x_hat(Nc,:),1,tf0(n,Nc),L,M,n);
%             A_hat_temp(indw)=Compute_A_DF2(x_hat(Nc,:),1,tf1,L,M,n);
%         end
%         A_hat_temp = sort(A_hat_temp);
%         A_hat(n,Nc) = A_hat_temp(3);
    end
end
A_hat = 2*A_hat; %% modified by DF (x2 only for real valued signals)



