clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 7 of the pper "Instantaneous Frequency and Amplitude
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
addpath(strcat([folder 'tools']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

% Amplitude modulated
amp0=((0.5.*cos(2*pi*(0:N-1)/(2*N/3)))+1)';

figure(1)
subplot(1,3,1)
colormap(flipud(gray))
plot(amp0)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 15, 'FontWeight', 'bold')
ylabel('Component amplitude','FontSize', 15, 'FontWeight', 'bold')


X(:,1) = amp0.*real(fmconst(N,0.3));
X(:,2) = amp0.*real(fmlin(N,0.15,0.3));


x0 = sum(X,2);

[tfr]  = tfrgab2(X(:,1), M, L);
Spect=(abs(tfr(1:round(M/2),:)));

figure(1)
subplot(1,3,2)
imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')
axis square





[tfr]  = tfrgab2(X(:,2), M, L);
Spect=(abs(tfr(1:round(M/2),:)));


figure(1)
subplot(1,3,3)
imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')
axis square


