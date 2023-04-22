clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative robustness evaluation for mask estimation only
%  using the RMSE to reproduce the results in Fig 2
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 16-feb-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'SSA']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'RD']));
addpath(strcat([folder '3DRD']));

snr_range = [-20 20]; % SNR range to compute
MCrep =10;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
% M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length


X(:,1) = real(fmlin(N,0.15,0.3));
X(:,2) = real(fmsin(N,0.3,0.45,700,1,0.3,+1));

x0 = sum(X,2);
Ncomp = size(X,2);                          %% number of components
X = real(transpose(X));

%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 5; % variance of the random walk in the temporal model
alpha = 0.5;
beta = 0.5;
detect = 0;

%% Calcul normalisation fft
A = 1/(sqrt(2*pi)*L);
C = -1 / (2*L^2);
K = 2 * L * sqrt(2*log(1/10^(-4)));  %% window length in samples
k = (-round(K/2):round(K/2)).^2;
g = A * exp( C * k);
G = sum(abs(fft(g)))/2;


%% Initialization
methods_name = {'Brevdo',...
                '3DRD [21]',...
                'PB',...
                'RD', ... 
                'Proposed EM-Lap',...
                'Proposed EM-TV'
                };
            
methods_to_use = [1 3 4 5 6];   % insert here the indices of the methods to compare (names above)
nb_methods = length(methods_to_use);
mm = [500 1000 2000];
ctime = zeros(MCrep,1);
Ctime = zeros(length(methods_to_use),length(mm));


%% Compute RMSE

for ind_met = 1:length(methods_to_use)

    for ind_m = 1:length(mm)
        M = mm(ind_m);

    for it = 1:MCrep   %% iterations
        clc;
        disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
        disp(strcat(['M : ', num2str(mm(ind_m))]))
        disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))


        % Add noise
        x = sigmerge(x0, randn(size(x0)), 10);
        [tfr]  = tfrgab2(x, M, L);
        Spect = abs(tfr(1:M/2,:)).^2;
        
        switch(methods_to_use(ind_met))
            case 1  %%Brevdo STFT
                tic
                [~,mask,tf] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                ctime(it) = toc;
            
            case 2  %% 3D RD
                tic
                   lambda_1 = 0.001; % in [0.001;0.1]
                   lambda_2 =  0.001; % in [0.001;0.1]
                   delta_1 = 2;
                   delta_2 = 2;
                   tf = ChirprateRD(x,M,Ncomp,lambda_1,lambda_2,delta_1,delta_2);
                ctime(it) = toc

            case 3  %% PB
                beta  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                div   = 2;   % 3 = Renyi
                PneiMask = 5;
                tic
                [mask,tf] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,0,PneiMask);
                ctime(it) = toc;

            case 4  %% RD
                 Nr = Ncomp;
                sigma_s = 0.09;
                clwin = 10;
                tic
                [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT, Cs_simple] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                ctime(it) = toc;
              
            case 5 %% EM LAP
                step_r = 10;
                step_v = 10;
                [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                tic
                [~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-3,step_r,step_v,ifplot,0,G);
                ctime(it) = toc;

             case 6 %% EM TV
                step_r = 10;
                step_v = 10;
                [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                tic
                [~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'TV',1e-4,step_r,step_v,ifplot,0,G);
                ctime(it) = toc;

        end  %% switch
    end %% repetitions


    Ctime(ind_met,ind_m) = mean(ctime);
    end
end  %% methods


for ind_met = 1:length(methods_to_use)
    for ind_m = 1:length(mm)
        disp(strcat([methods_name(ind_met),', M = ',num2str(mm(ind_m)),' : ',num2str(Ctime(ind_met,ind_m)) ]))
    end
end
 



