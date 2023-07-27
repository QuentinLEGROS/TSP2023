clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 2 of the pper "Instantaneous Frequency and Amplitude
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
addpath(strcat([folder 'Modul_EM']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

X(:,1) = real(fmlin(N,0.1,0.4));
X(:,2) = real(fmlin(N,0.3,0.15));


x0 = sum(X,2);

Ncomp = size(X,2);                          %% number of components

n_pad = 30;
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach


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
step_v = 4*step_r;


%% Compute ground truth
[tfr]  = tfrgab2(x0, M, L);
spect0=(abs(tfr(1:round(M/2),:)));

X = real(transpose(X));

% Add noise
x = sigmerge(x0, randn(size(x0)), 10);
[tfr]  = tfrgab2(x, M, L);
Spect = abs(tfr(1:M/2,:)).^2;
            

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-2,step_r,step_v,ifplot,0,0,0,0);
tf = min(max(1,tf-1),N/2);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
tf = tf(:,I);


figure(1)
subplot(2,2,2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-2,step_r,step_v,ifplot,1,0,0,0);
tf = min(max(1,tf-1),N/2);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
tf = tf(:,I);
            
figure(1)
subplot(2,2,1)
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
title('With interpolation')
axis square       
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear X;   
X(:,1) = real(fmconst(N,0.25));
X(:,2) = real(fmsin(N,0.2,0.3,700));

x0 = sum(X,2);
Ncomp = size(X,2);                          %% number of components

n_pad = 30;
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach

% Pour le mélange sinusoide + FM chirp
step_r = 4;
step_v = 10;


%% Compute ground truth
[tfr]  = tfrgab2(x0, M, L);
spect0=(abs(tfr(1:round(M/2),:)));

X = real(transpose(X));

% Add noise
x = sigmerge(x0, randn(size(x0)), 10);
[tfr]  = tfrgab2(x, M, L);
Spect = abs(tfr(1:M/2,:)).^2;
            

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-2,step_r,step_v,ifplot,0,0,0,0);
tf = min(max(1,tf-1),N/2);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
tf = tf(:,I);


figure(1)
subplot(2,2,4)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-2,step_r,step_v,ifplot,1,0,0,0);
tf = min(max(1,tf-1),N/2);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
tf = tf(:,I);
            
figure(1)
subplot(2,2,3)
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
title('With interpolation')
axis square       
            
            