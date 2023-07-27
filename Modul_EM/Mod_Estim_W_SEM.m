function [W_out,P,Mu_out,Amp_out]=Mod_Estim_W_SEM(Y,Fct,Ns,step_Nx,c,step_r,step_v,ifplot,bolpol,bolamp,G,N_samp,L)
    

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

% Y = Y./sum(Y,2);

%% Initialization
[N,M]=size(Y);% nb of time bins x nb of frequency bins
Nx=size(Fct,2); % nb of admissible frequencies
Nz = Nx/step_Nx; % subsampling the frequency grid for first iterations -> speed up the estimation
% ind0=cell(N,1); % cell for non-empty time bins
if ifplot
    cols = {'r-.', 'g-.', 'b-.', 'c--', 'm-.', 'g-x', 'w-o'};
end


wc=0.01*ones(N,Ns); % Initialization current W
muc=(floor(M/(Ns+1)):floor(M/(Ns+1)):M-floor(M/(Ns+1))).*ones(N,Ns); % Initialization current mu
    
iteEM = 0; % Number of EM iteration used to average W
pT=ones(1,Nz)./Nx; % Initialization uniform prior for mu (first EM iterations)

alpha=1.01*ones(1,Ns+1); % Initialization W priors fixed : alpha,beta (first EM iterations)


% Initialization without priors
C=Mod_comp_plik_multi(Y,Fct,wc,Ns); 
C=C-max(C,[],2)*ones(1,Nz);
P=exp(C);
P=P./(sum(P,2)*ones(1,Nz)); % Normalization
muc=denSampling(1:Nz,P); % Gibbs sampling
muc = muc.*step_Nx;
   
W_out=zeros(N,Ns); % Initialization output 
Mu_out=zeros(N,Ns); % Initialization output 
m_compt=1; % Initialization of EM iteration Count
Stopping_EM = 0; % EM stopping criterion
err0 = 10e9;
derr = 10e9;
errold = 10e9;

while ((m_compt<=20 && derr>1e-3 && err0>1e-10) || m_compt<=20)
    %% E-step
%     disp(['EM iteration : ',num2str(m_compt)])
    if m_compt<2 % first iterations without regularization
       C=Mod_comp_plik_multi(Y,Fct,wc,Ns); 
       C=C+ones(N,1)*log(pT);
       C=C-max(C,[],2)*ones(1,Nz);
       P=exp(C);
       P=P./(sum(P,2)*ones(1,Nz));
       muc = compute_itSMAP(P,Ns,step_r,step_v);
    else % Apply Lap prior
       C=Mod_comp_plik_multi(Y,Fct,wc,Ns);  % using SEM
       [P,muc]=compute_Sstep_TV(C,muc,c,N_samp,Ns,step_r,step_v);
    end

    
    %% Maximization step
    wc=Mod_Mstep_multi(Y,Fct,wc,alpha,Ns,muc);
    derr = norm(err0-errold);
    errold = err0;
    

    

    %% Plots
    if ifplot
        figure(2)
        subplot(3,1,1)
        imagesc(Y')
        hold on
        for ns = 1:Ns 
          IF = round(muc(:,ns));
          h(ns) = plot(1:N,IF, cols{ns});
          label{ns} = sprintf('mode %d', ns);
        end
        legend(h, label);
        xlabel('time [s]')
        ylabel('frequency [Hz]')
        title('Data (noisy)')
        subplot(3,1,2)
        plot(muc)
        ylim([0,250])
        
        subplot(3,1,3)
        plot(wc)

        pause(0.5)
    end
    
    Mu_out = muc;
    m_compt = m_compt+1;

    if (m_compt==10)  % Begin the last 5 iterations
        Stopping_EM = 1;
    end
    if Stopping_EM == 1
        W_out=W_out+wc;
        iteEM = iteEM + 1;
    end
end

%% Compute M
Mu_out = compute_itSMAP(P,Ns,step_r,step_v);

[~,nn]=sort(mean(Mu_out,1));
Mu_out = Mu_out(:,flip(nn));


W_out = W_out ./ iteEM; % mean of the (N_iter-N_bi) iterations

if bolpol
    Mu_out = poly_int(muc,step_r);
%     Mu_out = poly_in(muc);
end
Mu_out = min(max(Mu_out,1),Nx);

if ifplot
    %close all
    figure(2)
    imagesc(Y')
    hold on
    for c = 1:Ns 
      IF = round(muc(:,c));
      h(c) = plot(1:N,IF, cols{c});
      label{c} = sprintf('mode %d', c);
    end
    legend(h, label);
    xlabel('time [s]')
    ylabel('frequency [Hz]')
    title('Final estimation')
    pause(0.01)
end



if bolamp
    Amp_out = Estim_Amp(Mu_out,Y,W_out,G,M,L);
    Amp_out = sqrt(Amp_out);
end








