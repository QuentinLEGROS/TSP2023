clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 3 of the pper "Instantaneous Frequency and Amplitude
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
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

X(:,1) = real(fmlin(N,0.15,0.3));
X(:,2) = real(fmsin(N,0.3,0.45,700,1,0.3,+1));

x0 = sum(X,2);
tfr  = tfrgab2(x0, M, L);
Spect=(abs(tfr(1:round(M/2),:)));


figure
colormap(flipud(gray))
imagesc(Spect)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 12, 'FontWeight', 'bold')
ylabel('Normalized frequency','FontSize', 12, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})
