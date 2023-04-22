clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative robustness evaluation for mask estimation only
%  using the RMSE to reproduce the results in Fig 2
%  
%
%  Authors : Q.Legros, S. Meignen and D.Fourer
%  Date    : 01-feb-2023
%

%% Time-frequency representation parameters
N     = 500;         %% signal length
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
X(:,1) = (fmlin(N,0.1,0.4));
X(:,2) = (fmlin(N,0.3,0.15));

x0 = sum(X,2);
Ncomp = size(X,2);                          %% number of components

% Add noise
SNRi = inf;
x = sigmerge(x0, randn(size(x0)), SNRi);


%% Method parameters
lambda_1 = 0.001; % in [0.001;0.1]
lambda_2 =  0.001; % in [0.001;0.1]
delta_1 = 2;
delta_2 = 2;

%% Main algorithm
tf = ChirprateRD(x,M,Ncomp,lambda_1,lambda_2,delta_1,delta_2);


%% Plots
[tfr]  = tfrgab2(x, M, L);
Spect = abs(tfr(1:M/2,:)).^2;
figure(2)
imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')
hold on
plot(tf)
hold off
title('Without interpolation')
axis square







