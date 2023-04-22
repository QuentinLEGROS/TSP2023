clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Figure 10 of the pper "Instantaneous Frequency and Amplitude
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
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'RD']));
addpath(strcat([folder 'mData']));

addpath('./RRP_alg/');
addpath('./test/');
addpath('./TF_Toolbox/');

% Tfr parameter
L = 40;
M = 512;
M2 = floor(M/2);

[s,Fs] = audioread('/mData/OBVI.wav'); 

Fs_out = 5000;
s = resample(s,Fs_out,Fs);
N = length(s);

[tfr,stfr]  = tfrsgab2(s, M, L);
Spect = abs(tfr(1:M2,:)).^2;

Spectplot = abs(tfr(1:M2,:));

N = size(Spect,2);
t = (0:N-1)/Fs_out;
f = m_axis(M)/M*Fs_out;
f = f(1:M/2);

%% Method parameters

Ncomp = 2;                                 % Number of components
ifplot = 0;
Pnei = 8;
step_r = 5;
step_v = 1;

[Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-1,step_r,step_v,ifplot,0,0);

%% Plot
cols = {'r-', 'g-', 'b-', 'c-', 'm-', 'k-', 'w-'};

figure(1)
imagesc(t,f(1:M2),Spectplot);
% imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time (s)','FontSize', 15, 'FontWeight', 'bold')
ylabel({'Normalized';'Frequency'},'FontSize', 15, 'FontWeight', 'bold')
ylim([0 1000])
xlim([0.05,0.4])

hold on
for c = 1:Ncomp
  h(c) = plot(t,tf(:,c)/M*Fs_out, cols{c});
  label{c} = sprintf('mode %d', c);
end
yticks([1 500 1000])
yticklabels({0,0.25,0.5})
legend(h, label);
axis square




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % RD
 
 
figure(2)
imagesc(t,f(1:M2),Spectplot);
% imagesc(Spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time (s)','FontSize', 15, 'FontWeight', 'bold')
ylabel({'Normalized';'Frequency'},'FontSize', 15, 'FontWeight', 'bold')
ylim([0 1000])
xlim([0.05,0.4])


sigma_s = 0.09;
clwin = 10;
Nr = 2;
len_s = length(s);
x_hilbert = hilbert(s);
% x_hilbert = add_noise(x_hilbert, -2);
sigma_s = 0.07;
[STFT, TFR] = sst2(x_hilbert, sigma_s, M);
P = 1 - 1E-6;
[R_out, ~] = RRP_RD(STFT, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 'Nr', Nr);
tf = zeros(2164, 2);
for m_spl=1:Nr
    tf(:, m_spl) = fnval(R_out(m_spl).spline, (0:len_s-1)/len_s);
end
tf = tf*M/len_s; % normalization



cols = {'r-', 'g-', 'b-', 'k-', 'm-x', 'g-x', 'w-o'};
figure(2)
hold on
for c = 1:Ncomp
  h(c) = plot(t,tf(:,c)/M*Fs_out, cols{c});
  label{c} = sprintf('mode %d', c);
end
yticks([1 500 1000])
yticklabels({0,0.25,0.5})
legend(h, label);
axis square

