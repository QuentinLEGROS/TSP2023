function [Amp_out]=Oracle_EM(Y,Fct,Ns,G,tf0,M,L)
    

% Main algorithm: estimate the mixture weights
% 
% INPUT:
% Y         : Spectrogram of the MCS
% Fct       : Convolution matrix of the Gaussian kernel - data distribution
% Ns        : Number of spectral component (1 since mono-componant case)
% step_Nx   : Likelihood subsampling factor (for future research, 1 for now)
% reg_mu    : Spatial regularisation (between sequential time slices)
% c         : TV spatial prior weight
% cl        : Laplacian spatial prior weight
%
% OUTPUT:
% W_out    : Mixture weight estimates
% P        : IF Posterior distribution
% ind0     : Non zeros TF instants
%
% Author: Q.Legros

%% Initialization
[N,~]=size(Y);% nb of time bins x nb of frequency bins

wc=0.01*ones(N,Ns); % Initialization current W
iteEM = 0; % Number of EM iteration used to average W
alpha=1.01*ones(1,Ns+1); % Initialization W priors fixed : alpha,beta (first EM iterations)


W_out=zeros(N,Ns); % Initialization output 
m_compt=1; % Initialization of EM iteration Count
Stopping_EM = 0; % EM stopping criterion


% while (m_compt<=20)


    %% Maximization step
        W_out=Mod_Mstep_multi(Y,Fct,wc,alpha,Ns,tf0);
% 
%     
% 
%     m_compt = m_compt+1;
% 
%     if (m_compt==15)  % Begin the last 5 iterations
%         Stopping_EM = 1;
%     end
%     if Stopping_EM == 1
%         W_out=W_out+wc;
%         iteEM = iteEM + 1;
%     end
%     
%     figure(1)
%     plot(wc)
%     pause(0.5)
%     
%     
% end
% W_out = W_out ./ iteEM; 

Amp_out = Estim_Amp(tf0,Y,W_out,G,M,L);
Amp_out = sqrt(Amp_out);




