% Example 2 - bat signal analysis
clear all
close all
clc

folder = './';




%% required paths 
addpath(folder);
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'RD']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));

addpath('./RRP_alg/');
addpath('./test/');
addpath('./TF_Toolbox/');

load('bat2.mat');

Ncomp = 3;
N = length(x);
L = 9;
M = 512;
M2 = floor(M/2);
t = ((1:N)-1)/Fs*1000; %converted in ms
f = m_axis(M)/M*Fs;


[tfr,stfr]  = tfrsgab2(x, M, L);
spect=abs(tfr(1:M2,:)).^1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Lap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
imagesc(t,f(1:M2),spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')


ifplot = 0;
Pnei = 8;
Spect = abs(tfr(1:M/2,:)).^2;
step_r = 30;
step_v = 1;


[Fc]=comp_Fc(M,L);
Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-1,step_r,step_v,ifplot,1,0);
[mask] = compMask(round(tf),Pnei,N,0);

cols = {'r-', 'g-', 'b-', 'k-', 'm-x', 'g-x', 'w-o'};
figure(1)
hold on
for c = 1:Ncomp
  [ IF ] = mask2if( mask(1:M2,:,c) );
  h(c) = plot(t,IF/M*Fs, cols{c});
  
  label{c} = sprintf('mode %d', c);
    
end
legend(h, label);
title('Quentin (EM)')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % RD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
figure(2)
imagesc(t,f(1:M2)./M*Fs,spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 15, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 15, 'FontWeight', 'bold')

 Nr = Ncomp;
sigma_s = 0.09;
clwin = 10;
x_hilbert = hilbert(x);
sigma_s = 0.07;
[STFT, TFR] = sst2(x_hilbert, sigma_s, M);
P = 1 - 1E-6;
[R_out, ~] = RRP_RD(STFT, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 'Nr', Nr);
len_x = length(x);
tf = zeros(400, 3);
for m_spl=1:Nr
    tf(:, m_spl) = fnval(R_out(m_spl).spline, (0:len_x-1)/len_x);
end
tf = tf*M/len_x; % normalization
tf = tf*f(2); % mise a l'echelle

cols = {'r-', 'g-', 'b-', 'k-', 'm-x', 'g-x', 'w-o'};
figure(2)
hold on
for c = 1:Ncomp
  h(c) = plot(t,tf(:,c)/M*Fs, cols{c});
  label{c} = sprintf('mode %d', c);
    
end
yticks([1 f(M2/2)./M*Fs f(M2)./M*Fs])
yticklabels({0,0.25,0.5})
legend(h, label);
title('RD')
axis square