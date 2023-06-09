close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

phi1 = 1400*t+1350/(2*pi)*cos(2*pi*t + pi/2);
phi2 = 3400*t+550/(2*pi)*cos(2*pi*t);
IF1 = 1400-1350*sin(2*pi*t + pi/2);
IF2 = 3400-550*sin(2*pi*t);
s1_cos_lim = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_cos_lim = transpose([s1_cos_lim, s2]);
IFs_cos_lim = transpose([IF1, IF2]);
sigma_cos_lim = 0.0175;

Nr = 2;
Nfft = 512;
cas = 1;
poly_degree = 5;

%% SNR (RD + MR) computation
clwin = 10;
SNR_IN = -10:2:4;
NRep = 30;

[SNR_NEW_LCR, SNR_NEW_MR, SNR_S_MR, SNR_MB_MR, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD] =...
    test_RD_MR(modes_cos_lim, IFs_cos_lim, clwin, sigma_cos_lim, Nfft, poly_degree, SNR_IN, NRep);

plot_SNR_modes(SNR_IN, SNR_NEW_LCR, SNR_NEW_MR, SNR_S_MR, SNR_MB_MR);
savefig('F5_MR_cos');
saveas(gcf,'F5_MR_cos','epsc');
% close all;

plot_SNR_IFs(SNR_IN, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD);
savefig('F4_RD_cos');
saveas(gcf,'F4_RD_cos','epsc');
% close all;
