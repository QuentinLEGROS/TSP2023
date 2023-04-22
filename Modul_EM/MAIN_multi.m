clear all                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     % clear all
close all
clc
warning('off','all')

folder = './';
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'synchrosqueezedSTFT']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length
% X(:,1)    = (fmconst(N, 0.15));
% X(:,2)    = (fmconst(N, 0.35));
X(:,1) = fmlin(N,0.1,0.4);
X(:,2) = fmlin(N,0.4,0.1);
% X(:,1) = fmsin(N,0.3,0.45,320,1,0.3,+1);

x0 = sum(X,2);
Ncomp = size(X,2);                  %% number of components

SNR = inf;
% add noise
x = sigmerge(x0, randn(size(x0)), SNR);

%% Ground truth
[tfr]  = tfrgab2(x0, M, L);
Y0 = abs(tfr(1:M/2,:)).^2;

%% Data
[tfr]  = tfrgab2(x, M, L);
Y = abs(tfr(1:M/2,:)).^2;


[Fc]=comp_Fc(M,L);                  %% Data distribution
Fc = Fc + eps;
% Fc = repmat(Fc,[1,1,Ncomp]);

%% Initialize Parameters
c=1e-4; % TV regularization parameter
cl=1e-2; % Laplacian spatial prior reg hyperparameter
step_Nx = 1; % Depth grid subsampling
stepg = 1e-3; % Gradient step for the M-step
ifplot = 1;

%% Prior choice
% reg_mu='None';
reg_mu='Lap';
% reg_mu='TV';

%% EM algorithm
tic
% [a_out,T_out]=EM_ranging(Y',Fc,Ncomp,c,si,step_Nx,stepg,nb_label);
[W_out,P,tf]=Mod_Estim_W_EM_multi(Y',Fc,Ncomp,step_Nx,reg_mu,c,cl,ifplot);
toc
%% Plot reconstructions

%     figure(1)
% hold on
% plot(1:500,tf(:,1),'--')
%     plot(1:500,tf(:,2),'--')
