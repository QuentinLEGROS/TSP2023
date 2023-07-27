clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 1 of the pper "Instantaneous Frequency and Amplitude
%  Estimation in Multi-Component Signals Using Stochastic EM Algorithm"
%  
%
%  Authors : Q.Legros (legros.quentin2@hotmail.fr), D.Fourer, S.Meignen and
%  M.A.Colominas
%  Date    : 22-Avr-2023
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'PseudoBay']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

X(1:250,1) = real(fmconst(250,0.25));
X(251:N,1) = real(fmconst(250,0.3));

x0 = sum(X,2);

Ncomp = size(X,2);                          %% number of components

n_pad = 30;
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
step_r = 3;
step_v = 10;


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

% Add noise
SNR = 10;
x = sigmerge(x0, randn(size(x0)),SNR);
[tfr]  = tfrgab2(x, M, L);
Spect = abs(tfr(1:M/2,:)).^2;
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'TV',1e1,step_r,step_v,ifplot,0,0,0,0);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
x_hat = x_hat(I,:);
tf = tf(:,I);



figure(1)
subplot(1,2,1)
imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 15, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 15, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})
hold on
plot(tf)
hold off
title('SNR = 10dB')
axis square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add noise
SNR = -5;
x = sigmerge(x0, randn(size(x0)),SNR);
[tfr]  = tfrgab2(x, M, L);
Spect = abs(tfr(1:M/2,:)).^2;

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'TV',1e1,step_r,step_v,ifplot,0,0,0,0);
[mask] = compMask(round(tf),Pnei,N,0);
x_hat = zeros(Ncomp,N);
for c = 1:Ncomp
   x_hat(c,:) = 2*real(rectfrgab(tfr .* mask(1:M,:,c), L, M)); 
end

[I,~] = match_components(X, x_hat); 
x_hat = x_hat(I,:);
tf = tf(:,I);
            
figure(1)
subplot(1,2,2)
imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 15, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 15, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})
hold on
plot(tf)
hold off
title('SNR = -5dB')
axis square       