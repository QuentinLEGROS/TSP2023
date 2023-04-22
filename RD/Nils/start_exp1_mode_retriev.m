clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative evaluaton for mask estimation and 
%  mode retrieval considering a multicomponent signal made 
%  of 3 components merged with a white Gaussian noise
%  
%
%  Authors : D. Fourer (dominique@fourer.fr) and Q. Legros 
%  Date    : 13-feb-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'SSA']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

if ~exist('./figs', 'dir')
  mkdir('figs');    
end

%%%       generate signals, oracle
snr_in = 15;


N=500;
Ncomp = 3;

X = zeros(N,Ncomp);
X(:,1) = real(fmconst(N, 0.1));
X(:,2) = real(fmlin(N,0.13,0.3));
X(:,3) = real(fmsin(N,0.3,0.45,320,1,0.3,+1));
x0 = sum(X,2);


% TFR - parameters
M       = 500;       %% nombre de bins frequentiels
mm      = m_axis(M);
L       = 20;        %% taille de la fenetre d'analyse en samples

gamma_mask = 0.12; % 0.1;  %% threshold for the mask
mask_ref = zeros(M,N, Ncomp);
for nc = 1:Ncomp
  [tfr_tmp] = tfrgab2(X(:,nc), M, L);
  mask_ref(:,:, nc) = oracle_mask(tfr_tmp, gamma_mask);
  %imagesc(mask_ref(:,:, nc))
end


methods_name = {'Brevdo SST', 'Brevdo STFT', 'SSA', 'Proposed \beta=1, Beta div.', 'Proposed KL div.', ... %% 1-5
                'Proposed \alpha=0.9, Renyi div.', sprintf('Oracle \\Gamma=%.2f', gamma_mask) ...          %% 6-7
                'Nils 1', 'Nils 2', 'Nils 3', 'Nils 4', 'Nils spl', 'Nils spl lcr'};                       %% 8-13

%methods_to_use = [1 4 6  8 9 12 13];   % insert here the indices of the methods to compare (names above)
methods_to_use = [1  8 9 12 13];   % insert here the indices of the methods to compare (names above)


nb_methods = length(methods_to_use);


x = sigmerge(x0, randn(size(x0)), snr_in); %% plus precis  
x = x - mean(x);

%% step 1, display TFR (STFT and SST)
[tfr,stfr] = tfrsgab2(x, M, L);

figure(111)
m_max = round(M/2);
f = ((1:m_max)-1)/M;
t = ((1:length(x))-1);
imagesc(t,f, abs(tfr(1:m_max,:)).^0.4);
set(gca,'Ydir','normal');
colormap gray;colormap(flipud(colormap));
xlabel('time index', 'FontSize', 15);
ylabel('normalized frequency', 'FontSize', 15);
title(sprintf('spectrogram, SNR=%.2f dB', snr_in), 'FontSize', 15);
saveas(gcf, sprintf('figs/spectrogram.eps'), 'epsc'); 

figure(112)
m_max = round(M/2);
f = ((1:m_max)-1)/M;
t = ((1:length(x))-1);
imagesc(t,f, abs(stfr(1:m_max,:)).^0.4);
set(gca,'Ydir','normal');
colormap gray;colormap(flipud(colormap));
xlabel('time index', 'FontSize', 15);
ylabel('normalized frequency', 'FontSize', 15);
title(sprintf('squared modulus of the synchrosqueezed STFT, SNR=%.2f dB', snr_in), 'FontSize', 15);
saveas(gcf, sprintf('figs/sst.eps'), 'epsc'); 


for ind_met = 1:nb_methods
        
    fprintf(1,'Computing Method [%d]: %s ...\n', methods_to_use(ind_met), methods_name{methods_to_use(ind_met)});
	RQF_out = zeros(Ncomp, nb_methods);
              
    switch(methods_to_use(ind_met))
    	case 1  %%Brevdo SST
            [tfr,stfr] = tfrsgab2(x, M, L);   %% compute SST
            %% K=5
            K = 8; %8; %%vicinity of the ridge
            [~, mask] = Brevdo_modeExtract(stfr, L, Ncomp, K);
            
            x_hat = zeros(N,Ncomp);
            for c = 1:Ncomp
               x_hat(:,c) = real(rectfrsgab(stfr .* mask(:,:,c), L, M));       
            end
                                 
        case 2  %%Brevdo STFT
            [tfr] = tfrgab2(x, M, L);   %% compute SST
            K = 8; %8; %%vicinity of the ridge
            [~, mask] = Brevdo_modeExtract(tfr, L, Ncomp, K);
                
            x_hat = zeros(N,Ncomp);
            for c = 1:Ncomp
               x_hat(:,c) = real(rectfrgab(tfr .* mask(:,:,c), L, M));       
            end
        case 3  %% SSA
            nc = Ncomp+1;             %% number of components
            epsilon = 1e-2;     %% singular spectrum thresholding parameter (increase the value to improve the robustness to noise)

            x_tmp = ssa_decomp(x, L, nc, epsilon);
            x_hat = x_tmp(:,1:Ncomp);
            
            
            mask = zeros(M,N,Ncomp);
            for nc = 1:Ncomp
              [tfr_tmp] = tfrgab2(x_hat(:,nc), M, L);
              mask(:,:, nc) = oracle_mask(tfr_tmp, gamma_mask);       
            end

  %% Bayesian
                    %disp(strcat(['Method :', methods_name{3}]));
                        % Values to try for beta : [0.2-0.5-0.7-0.9] %  Should give better results
                        %                                               with lower beta when the SNR 
                        %                                               is low
                        % Values to try for alpha : [0.2-0.5-0.7-0.9] % I think 0.5 should be good,
                        %                                               lets see
                        %
                        % Do not touch to ds, It gave good results for ds=3, and I do not think
                        % this should change a lot of thing on the signals we will try
                        % 
                        % There is still something to check if you have the motivation : the
                        % variance of the Gaussian used for the ridge removal. As discussed, that
                        % might be normal and just the result of a specific behavious that we
                        % haven't understood yet. However, it would be very cool to be able to reproduce 
                        % the Gaussian corresponding to the signal at a given time without any
                        % noise and for only one ridge. I will also check when I will have some
                        % time. Not before monday for me.
                        
        case 4  %% Beta divergence
          %% Bayesian method parameters
          beta  = 1;   % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
          alpha = 0.9; % Renyi divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 1
          ds    = 2;    % variance of the random walk in the temporal model
          div   = 2;   % 1 = KL
                                     % 2 = beta
                                     % 3 = Renyi
          ifplot =  0;
          Pnei   =  10;
          [tfr,stfr]  = tfrsgab2(x, M, L);
          % pseudoBay(tfr,Ncomp,M,L,div,beta,alpha,ds,Pnei,ifplot)
          %[mask] = pseudoBay(tfr, Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
          [mask, rectfr] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);              
          x_hat = zeros(N,Ncomp);
          for c = 1:Ncomp
            %x_hat(:,c) = real(rectfrsgab(tfr .* mask(:,:,c), L, M));
            x_hat(:,c) = real(rectfrgab(rectfr(:,:,c), L, M));
            %x_hat(:,c) = 2*real(rectfrgab(tfr(1:Mh,:) .* mask(1:Mh,:,c), L, M, mm(1:Mh))).';
            %% SST (uncomment for SST version)
            %x_hat(:,c) = real(rectfrsgab(stfr .* mask(:,:,c), L, M)).'; 
          end
                        
        case  5  %% KL divergence
          beta  = 1;   % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
          alpha = 0.9; % Renyi divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 1
          ds    = 3;    % variance of the random walk in the temporal model
          div   = 1;   % 1 = KL
                                   % 2 = beta
                                   % 3 = Renyi
          ifplot =  0;
          Pnei   =  10;
          [tfr,stfr]  = tfrsgab2(x, M, L);
          % pseudoBay(tfr,Ncomp,M,L,div,beta,alpha,ds,Pnei,ifplot)
          [mask, rectfr] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
          x_hat = zeros(N,Ncomp);
          for c = 1:Ncomp
            %x_hat(:,c) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
            x_hat(:,c) = real(rectfrgab(rectfr(:,:,c), L, M));
            %x_hat(:,c) = real(rectfrgab(tfr(1:Mh,:) .* mask(1:Mh,:,c), L, M, mm(1:Mh))).';
            %% SST (uncomment for SST version)
            %x_hat(:,c) = real(rectfrsgab(stfr .* mask(:,:,c), L, M)).'; 
          end
                    
        case 6  %% Renyi diverg (Bug)
          beta  = 1;   % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
          alpha = 0.9; % Renyi divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 1
          ds    = 3;    % variance of the random walk in the temporal model
          div   = 3;   % 1 = KL
                                   % 2 = beta
                                   % 3 = Renyi
          ifplot =  0;
          Pnei   =  10;
          [tfr,stfr]  = tfrsgab2(x, M, L);
          % pseudoBay(tfr,Ncomp,M,L,div,beta,alpha,ds,Pnei,ifplot)
          %[mask] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
          [mask, rectfr] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
          x_hat = zeros(N,Ncomp);
          for c = 1:Ncomp
            %x_hat(:,c) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
            x_hat(:,c) = real(rectfrgab(rectfr(:,:,c), L, M));
            %x_hat(:,c) = 2*real(rectfrgab(tfr(1:Mh,:) .* mask(1:Mh,:,c), L, M, mm(1:Mh))).';
            %% SST (uncomment for SST version)
            %x_hat(:,c) = real(rectfrsgab(stfr .* mask(:,:,c), L, M)).'; 
          end
          
        case 7  %% ORACLE
          mask = mask_ref;
          x_hat = zeros(N,Ncomp);
          [tfr]  = tfrgab2(x, M, L);
          for c = 1:Ncomp
            x_hat(:,c) = real(rectfrgab(tfr .* mask_ref(:,:,c), L, M));
          end
          
        case 8  %% simple method
           sigma_s = 0.09;
           clwin = 10;
           [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M,Ncomp, sigma_s, clwin );
           for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_SR_Cl(c,:));
           end
           
        case 9 %% simple method2
            sigma_s = 0.09;
            clwin = 10;
            [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Ncomp, sigma_s, clwin );
            for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_SR_MB(c,:));
            end 
                    
        case 10 %% Linear chirp reconstruction method
            sigma_s = 0.09;
            clwin = 10;
            [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Ncomp, sigma_s, clwin );
            for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_LCR_Cl(c,:));
            end
            
        case 11 %% Linear chirp reconstruction method 2
            sigma_s = 0.09;
            clwin = 10;
            [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Ncomp, sigma_s, clwin );
            for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_LCR_MB(c,:));
            end
            
        case 12 %% Nils method
            sigma_s = 0.09;
            clwin = 10;
            smooth_p = 1-1e-10;
            [ m_SR, m_LCR, ~] = Nils_modeExtract2( x, M, Ncomp, sigma_s, clwin, smooth_p );
            for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_SR(c,:));
            end
             
        case 13 %% Nils method 2 LCR
            sigma_s = 0.09;
            clwin = 10;
            smooth_p = 1-1e-10;
            [ m_SR, m_LCR, ~] = Nils_modeExtract2( x, M, Ncomp, sigma_s, clwin, smooth_p );
            for c = 1:Ncomp
             x_hat(:,c) = 2*real(m_LCR(c,:));
            end
 
        otherwise 
          error('Invalid method'); 
    end  %% switch
  
  figure(methods_to_use(ind_met))
  plot_result(X.', x_hat.', methods_name{methods_to_use(ind_met)}, [-2 2]);  
  saveas(gcf, sprintf('figs/exp1_method%d.eps', methods_to_use(ind_met)), 'epsc');
  
  fprintf(1, '\n\n------------Method %s---------\n', methods_name{methods_to_use(ind_met)});
  [ I, ~ ] = match_components(X.', x_hat.');
  for nc = 1:Ncomp
   fprintf(1,'component %d :\t\t RQF=%.2f \t DICE=%.2f \n', nc, RQF(X(:,nc), x_hat(:,I(nc))), dice(mask_ref(:,:,nc),mask(:,:,I(nc))));   
  end
  
  
end

eps2pdf('./figs');







