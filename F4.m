clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 4 of the pper "Instantaneous Frequency and Amplitude
%  Estimation in Multi-Component Signals Using Stochastic EM Algorithm"
%  
%
%  Authors : Q.Legros (legros.quentin2@hotmail.fr), D.Fourer, S.Meignen and
%  M.A.Colominas
%  Date    : 22-Avr-2023
%
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'RD']));

snr_range = [-20 20]; % SNR range to compute
MCrep =50;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

X(:,1) = real(fmlin(N,0.15,0.3));
X(:,2) = real(fmsin(N,0.3,0.45,700,1,0.3,+1));


x0 = sum(X,2);

Ncomp = size(X,2);                          %% number of components

n_pad = 30;
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 5; % variance of the random walk in the temporal model
alpha = 0.5;
beta = 0.5;
detect = 0;
%% Compute ground truth
tf0=zeros(N,1);
for ns = 1:Ncomp
    [tfr]  = tfrgab2(X(:,ns), M, L);
    spect=(abs(tfr(1:round(M/2),:)));
    for i=1:N
        [~,mm]=max(spect(:,i));
        tf0(i,ns)=mm;
    end
end
X = real(transpose(X));


%% Calcul normalisation fft
A = 1/(sqrt(2*pi)*L);
C = -1 / (2*L^2);
K = 2 * L * sqrt(2*log(1/10^(-4)));  %% window length in samples
k = (-round(K/2):round(K/2)).^2;
g = A * exp( C * k);
G = sum(abs(fft(g)))/2;


%% Initialization
methods_name = {'Carmona [17]',...
                'PB-\alpha=0.5 10]',...
                'PB-\beta=0.5 [10]',...
                'RD [14]', ... 
                'Proposed EM-Lap',...
                'Proposed EM-TV'
                };
            
methods_to_use = [1 2 3 4 5 6];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):4:snr_range(2);

L2ErPos_out = zeros(length(SNRt), nb_methods);

%% Compute RMSE
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        L2ErPos_tmp = zeros(1,length(MCrep));

        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            

            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            [tfr]  = tfrgab2(x, M, L);
            Spect = abs(tfr(1:M/2,:)).^2;
            
            switch(methods_to_use(ind_met))
                case 1  %%Brevdo STFT
                    [~,mask,tf] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                    tf = tf';
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                       x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));       
                    end
                case 2  %% Alpha divergence
                    alpha  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 3;   % 3 = Renyi
                    [mask,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr(:,:) .* mask(:,:,c), L, M)); 
                    end

                case 3  %% Beta divergence
                    beta  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 2;   % 3 = Renyi
                    [mask,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr(:,:) .* mask(:,:,c), L, M)); 
                    end
                case 4  %% simple
                        Nr = Ncomp;
                        sigma_s = 0.09;
                        clwin = 10;
                        [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT, Cs_simple] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                        tf = Cs_simple';
                      [mask] = compMask(tf,Pnei,N,0);
                      x_hat = zeros(Ncomp,N);
                      for c = 1:Ncomp
                         x_hat(c,:) = real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
                      end
                  
                case 5 %% EM LAP
                    step_r = 10;
                    step_v = 10;
                    [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                    [~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-3,step_r,step_v,ifplot,0,0,G,1);
                    [mask] = compMask(tf,Pnei,N,0);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                       x_hat(c,:) = real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
                    end

                 case 6 %% EM TV
                    step_r = 10;
                    step_v = 10;
                    [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                    [~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'TV',1e-4,step_r,step_v,ifplot,0,0,G,1);
                    [mask] = compMask(tf,Pnei,N,0);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                       x_hat(c,:) = real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
                    end

            end  %% switch

             % Match components and reordering for comparison
            [I,~] = match_components(X, x_hat); 
            x_hat = x_hat(I,:);
            tf = tf(:,I);

            L2ErPos_tmp(it) = 10*log10(sum(sum((tf(n_pad:(end-n_pad),:)-tf0(n_pad:(end-n_pad),:)).^2))./(N*M*M));
        end    %% for methods
        L2ErPos_out(indsnr, ind_met) = mean(L2ErPos_tmp);

    end  %% methods
end %% snrs



%% Plot
cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-o' 'g-o' 'r-o' 'k--' 'b--' 'g--' 'r--' 'k-v' 'b-v' 'g-v'};
leg = {};
figure(1)
axis square
for ind_met =  1:nb_methods
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, L2ErPos_out(:,ind_met), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 15, 'bold')
 ylabel('RMSE (dB)', 'FontSize', 15)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

grid
