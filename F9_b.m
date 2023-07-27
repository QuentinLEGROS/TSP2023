clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 9b of the pper "Instantaneous Frequency and Amplitude
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
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'Compute_Amplitude_DF']));

snr_range = [-20 20]; % SNR range to compute
MCrep = 5;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

% Amplitude modulated
amp0(:,1)=(1:-0.5/(N-1):0.5);
amp0(:,2)=(0.5:0.5/(N-1):1);
% amp0(:,2)=((0.5.*cos(2*pi*(0:N-1)/(2*N/3)))+1)'

X(:,1) = amp0(:,1).*real(fmlin(N,0.1,0.4));
X(:,2) = amp0(:,2).*real(fmlin(N,0.3,0.15));


Ncomp = size(X,2);
x0 = sum(X,2);

%% Method parameters
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 10; % variance of the random walk in the temporal model
n_pad = 50; 
detect=0;
div = 4;
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 5;
alpha = 0.5;
beta = 0.5;


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

%% Calcul step_r and step_v
sigm_d = round(sqrt((M/2)/(pi*L)));
step_r = 3*sigm_d;
step_v = 3*step_r;

%% Initialization
methods_name = {'EM',...
                'Oracle EM',...
                'PB',...
                'Oracle PB',...
                'Local',...
                'Oracle Local',...
                };
            
            
methods_to_use = 2%[1 2 3 4 5 6];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):4:snr_range(2);

MAE_out = zeros(length(SNRt), nb_methods);

%% Compute RQF
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        MAE_tmp = zeros(MCrep,1);
        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
            
            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            switch(methods_to_use(ind_met))
                case 1  %% Proposed EM 
                        [tfr]  = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                        ifplot = 0;
                        [~,~,tf,Amp]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e0,step_r,step_v,ifplot,1,1,G,1,L);
                        tf = min(max(tf,1),M/2);
                        tf1(:,:,it) = tf;
                case 2  %% Oracle EM
                        [tfr]  = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                        [Amp]=Oracle_EM(Spect',Fc,Ncomp,G,tf0,M,L);
                        tf = tf0;
                case 3  %% Beta divergence
                        alpha  = 0.4;
                        beta   = 0.2;
                        [tfr]  = tfrgab2(x, M, L);
                        [tf,Amp] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,0,PneiMask,tf0);
                case 4
                        alpha  = 0.4;
                        beta   = 0.2;
                        [tfr]  = tfrgab2(x, M, L);
                        [tf,Amp] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,1,PneiMask,tf0);
                case 5 %% Méthode Dominique
                        [tfr]  = tfrgab2(x, M, L);
                        tf = squeeze(tf1(:,:,it));
                        [Amp] = Oracle_Amp_DF(tfr,Ncomp,M,L,tf); 
                case 6  %% Méthode Dominique Oracle
                        [tfr]  = tfrgab2(x, M, L);
                        tf = tf0;
                        [Amp] = Oracle_Amp_DF(tfr,Ncomp,M,L,tf); 
                        
            end  %% switch
            tf = min(max(tf,1),M/2);
            [mask] = compMask(tf,Pnei,M/2,0);
            x_hat = zeros(Ncomp,N);
            for c = 1:Ncomp
               x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
            end
            
            % Match components and reordering for comparison
            [I,~] = match_components(X, x_hat); 
            Amp = Amp(:,I);
            
            % MAE
            MAE_tmp(it) = sum(sum(abs(Amp(n_pad:end-n_pad,:) - amp0(n_pad:end-n_pad,:))))/(N*Ncomp);
        end    %% for methods
        MAE_out(indsnr, ind_met) = mean(MAE_tmp(~isnan(MAE_tmp)));
    end  %% methods
 
end %% snrs

% Normalization
MAE_out = MAE_out ./(N*Ncomp);
% 
figure
subplot(2,1,1)
plot(amp0)
subplot(2,1,2)
plot(Amp)


%% Plot
cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-s' 'g-.' 'r-.' 'b-v'  'b--' 'b*'};
leg = {};

figure(1)
iii = 0;
for ind_met =  1:nb_methods
 if ind_met == 1
  hold off
 else
  hold on
 end

 h(ind_met) = plot(SNRt, squeeze(MAE_out(:,ind_met)), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 15)
 ylabel('RMAE', 'FontSize', 15)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

grid
axis square
set(gca, 'YScale', 'log')
