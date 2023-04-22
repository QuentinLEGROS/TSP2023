clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative robustness evaluaton for mask estimation and 
%  mode retrieval considering a signal made of 1 linear chirp merged 
%  with a white Gaussian noise
%  
%
%  Authors : D. Fourer (dominique@fourer.fr) and Q. Legros 
%  Date    : 13-feb-2021
%  updated : 5-jul-2021
%

folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'SSA']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'Nils']));



snr_range = [-20 20]; % SNR range to compute
MCrep = 100; %100;     % number of Monte Carlo iterations (increase to obtain smooth curves)


%% Load signal (linear chirp)
N     = 500;                        %% signal length
x0    = real(fmlin(N,0.13,0.3));    %% linear chirp
Ncomp = 1;                          %% number of components

n_pad = 15;  %% number of samples to ignore for RQF computation

% TFR - parameters
M       = 500;       %% nombre de bins frequentiels
L       = 20;        %% taille de la fenetre d'analyse en samples


gamma_mask = 0.12; % 0.1;  %% threshold for the mask
[tfr_ref, stfr_ref] = tfrsgab2(x0, M, L);
mask_ref            = oracle_mask(tfr_ref, gamma_mask);


methods_name = {'Brevdo SST', 'Brevdo STFT', 'SSA', 'Proposed \beta=1, Beta div.', 'Proposed KL div.', ... %% 1-5
                'Proposed \alpha=0.9, Renyi div.', sprintf('Oracle \\Gamma=%.2f', gamma_mask) ...          %% 6-7
                'Nils 1', 'Nils 2', 'Nils 3', 'Nils 4', 'Nils spl', 'Nils spl lcr'};                       %% 8-13

methods_to_use = [1 4 6  8 9 12 13];   % insert here the indices of the methods to compare (names above)


nb_methods = length(methods_to_use);

%SNRt  = sort(real(logspace(log10(snr_range(2)),log10(snr_range(1)), 50)),'ascend');
SNRt = snr_range(1):2:snr_range(2);


RQF_out  = zeros(length(SNRt), nb_methods);
Dice_out = zeros(length(SNRt), nb_methods);

for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);
  


    for ind_met = 1:length(methods_to_use)
        
        RQF_tmp  = zeros(1,length(MCrep));
        Dice_tmp = zeros(1,length(MCrep));
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
            
            % Add random noise noise
            %noise = randn(size(s))*std(s)/db2mag(SNRi);
            %x = s + noise;
            x = sigmerge(x0, randn(size(x0)), SNRi); %% plus precis  
                
                
            switch(methods_to_use(ind_met))
            	case 1  %%Brevdo SST
                     %disp(strcat(['Method :', methods_name{1}]));
                     [tfr,stfr] = tfrsgab2(x, M, L);   %% compute SST
                     %% K=5
                     K = 8; %8; %%vicinity of the ridge
                     [~, mask] = Brevdo_modeExtract(stfr, L, Ncomp, K);
                     [ x_hat ] = real(rectfrsgab(stfr .* mask, L, M));
                     
            	case 2  %%Brevdo STFT
                     %disp(strcat(['Method :', methods_name{1}]));
                     [tfr] = tfrgab2(x, M, L);   %% compute SST
                     
                     K = 8; %8; %%vicinity of the ridge
                     [~, mask] = Brevdo_modeExtract(tfr, L, Ncomp, K);
                     [ x_hat ] = real(rectfrgab(tfr .* mask, L, M));       
                     
                case 3  %% SSA
                     %disp(strcat(['Method :', methods_name{2}]));
                     nc = Ncomp+1;             %% number of components
                     epsilon = 1e-2;     %% singular spectrum thresholding parameter (increase the value to improve the robustness to noise)
                    %L  = 20;

                    x_tmp = ssa_decomp(x, 20, nc, epsilon);
                    x_hat = x_tmp(:,1).';
                    
                    [tfr_tmp] = tfrgab2(x_hat, M, L);
                    mask = oracle_mask(tfr_tmp, gamma_mask);       


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
                        [tfr,stfr]  = tfrsgab2(x, M, L);
                        % (tfr,Ncomp,M,L,div,beta,alpha,ds,Pnei,ifplot)
                        [mask] = pseudoBay(tfr,Ncomp, M,  L, div, beta, alpha, ds, 15, ifplot);
                        %[mask] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, ifplot);
                        [ x_hat ] = real(rectfrgab(tfr .* mask, L, M));
                        %% SST (uncomment for SST version)
                        %[ x_hat ] = real(rectfrsgab(stfr .* mask, L, M)); 
                        
                case 5  %% KL divergence
                        beta  = 1;   % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha = 0.9; % Renyi divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 1
                        ds    = 2;    % variance of the random walk in the temporal model
                        div   = 1;   % 1 = KL
                                   % 2 = beta
                                   % 3 = Renyi
                        ifplot =  0;
                        [tfr,stfr]  = tfrsgab2(x, M, L);
                        [mask] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, 15, ifplot);
                        [ x_hat ] = real(rectfrgab(tfr .* mask, L, M));
                        %% SST (uncomment for SST version)
                        %[ x_hat ] = real(rectfrsgab(stfr .* mask, L, M));    
                    
                    
                case 6  %% Renyi diverg (Bug)
                        beta  = 1;   % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha = 0.9; % Renyi divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 1
                        ds    = 3;    % variance of the random walk in the temporal model
                        div   = 3;   % 1 = KL
                                   % 2 = beta
                                   % 3 = Renyi
                        ifplot =  0;
                        [tfr,stfr]  = tfrsgab2(x, M, L);
                        [mask] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, 15, ifplot);
                        [ x_hat ] = real(rectfrgab(tfr .* mask, L, M));
                        %% SST (uncomment for SST version)
                        %[ x_hat ] = real(rectfrsgab(stfr .* mask, L, M)); 
                        
                case  7  %% ORACLE
                  mask = mask_ref;
                  [tfr]  = tfrgab2(x, M, L);
                  [ x_hat ] = real(rectfrgab(tfr .* mask, L, M));
                  
                case 8  %% simple method
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                  [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                  x_hat = 2 * real(m_SR_Cl);
        
                case 9 %% simple method2
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                  [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                  x_hat = 2 * real(m_SR_MB);  
                    
                case 10 %% Linear chirp reconstruction method
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                  [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                  x_hat = 2 * real(m_LCR_Cl);    
                    
                case 11 %% Linear chirp reconstruction method 2
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                  [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                  x_hat = 2 * real(m_LCR_MB); 

                case 12 %% Nils method
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                   smooth_p = 1-1e-5; %1-1e-10
                   [ m_SR, m_LCR, ~] = Nils_modeExtract2( x, M, Nr, sigma_s, clwin, smooth_p );
                   x_hat = 2 * real(m_SR);  
                  
                case 13 %% Nils method 2 LCR
                   Nr = 1;
                   sigma_s = 0.09;
                   clwin = 10;
                   smooth_p = 1-1e-5;
                   [ m_SR, m_LCR, ~] = Nils_modeExtract2( x, M, Nr, sigma_s, clwin, smooth_p );
                   x_hat = 2 * real(m_LCR);
                  
                otherwise
                   error('Invalid method');
            end  %% switch
            RQF_tmp(it) = RQF(x0(n_pad:(end-n_pad)), x_hat(n_pad:(end-n_pad)).');
            Dice_tmp(it) = dice(mask, mask_ref);
            
            disp(strcat(['RQF : ', num2str(RQF_tmp(it)),'Dice : ', num2str(Dice_tmp(it))]));
        end    %% for methods
        RQF_out(indsnr, ind_met) = mean(RQF_tmp);
        Dice_out(indsnr, ind_met) = mean(Dice_tmp);
        disp(strcat(['RQF : ', num2str(RQF_out(indsnr, ind_met))]));
    end  %% methods
 
end %% snrs
    

cols         = {'k-x'  'b-o' 'g-s' 'b-^' 'g-.' 'r-.' 'r-v'   'b-.' 'b--' 'b*' 'k-.' 'o-^' 'v-x'};

h   = [];
leg = {};
figure(11)  %% RQF as a function of SNR
for ind_met = 1:nb_methods
 fprintf(1,'Method %s. min:%.2f  max:%.2f  median:%.2f \n', methods_name{ind_met}, min(RQF_out(:, ind_met)),max(RQF_out(:, ind_met)), median(min(RQF_out(:, ind_met))));
    
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, RQF_out(:,ind_met), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 14)
 ylabel('RQF (dB)', 'FontSize', 14)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end

legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)
grid
if ~exist('figs', 'dir')
  mkdir('./figs');
end
saveas(gcf, sprintf('figs/exp2_compare_RQF.eps'), 'epsc');
eps2pdf('./figs');


h   = [];
leg = {};
figure(22)  %% RQF as a function of SNR
for ind_met = 1:nb_methods
 fprintf(1,'Method %s. min:%.2f  max:%.2f  median:%.2f \n', methods_name{ind_met}, min(Dice_out(:, ind_met)),max(Dice_out(:, ind_met)), median(min(Dice_out(:, ind_met))));
    
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, Dice_out(:,ind_met), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 14)
 ylabel('Dice', 'FontSize', 14)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end

legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)
grid
if ~exist('figs', 'dir')
  mkdir('./figs');
end
saveas(gcf, sprintf('figs/exp2_compare_dice.eps'), 'epsc');
eps2pdf('./figs');

