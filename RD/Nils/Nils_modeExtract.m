function [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin )
% [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin )
%
%  Nils mode extraction method
%
% input:
% x : signal
% M : number of frequency bins
% Nr: number of components to extract
% sigma_s : (optional, default 0.09) window width
% clwin : (optional, default =10) frequency clearing window
%
% output (estimated modes):
%
% m_SR_Cl  : Classical method
% m_SR_MB  : 
% m_LCR_Cl : 
% m_LCR_MB : 
% STFT     : STFT matrix
%
%


N = length(x);

if ~exist('sigma_s', 'var')
   sigma_s = 0.09;  
end

if ~exist('clwin', 'var')
   clwin = 10;  
end

if ~exist('create_gaussian_window')
    addpath('./Nils');
end


Nfft = M;

% Nr = 1;
% smooth_p = 1 - 10^(-4);
% t = (0:N-1)'/N;
% 
% [g, Lh] = create_gaussian_window(N, Nfft, sigma_LC);
% [STFT, omega, ~, QM, ~, tau] = FM_operators(x, L, Nfft, g, Lh, sigma_LC);
% [Spl_LC, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_LC, smooth_p);
% 
% Spl = Spl_LC;
% sigma_s = sigma_LC;
% IF1 = fnval(Spl(1).spline, t);
% DF1 = fnval(fnder(Spl(1).spline), t);
% R1 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF1.^2);

        
[g, Lh] = create_gaussian_window(N, Nfft, sigma_s);
[STFT, omega, ~, QM, ~, tau] = FM_operators(x, N, Nfft, g, Lh, sigma_s);
% QM : chirp rate

%  fprintf('Classic, ');
%% [8] extraction ridge simple (classique ) [R. Carmona,  W. Hwang,  and B. Torresani,  ?Characterization of signals  by the ridges  of their wavelet  transforms,?IEEETransactions on Signal Processing, vol. 45, no. 10, pp. 2586?2590, Oct 1997]
   [Cs_simple] = exridge_mult(STFT, Nr, 0, 0, clwin);
   % CS_simple : ridge
        
%         Spl_Cl = struct('spline', cell(1, Nr));
%         for m=1:Nr
%             Spl_Cl(m).spline = spline((0:L-1)/L, (Cs_simple(m, :) - 1)*L/Nfft);
%         end
%         [m_SR_Cl, m_LCR_Cl, IF_Cl] = R1_MR_and_LCR_spl(STFT, Spl_Cl, g, Lh, sigma_s, Nr, Nfft, L);

%% reconstruction (m_SR_Cl : simple reconstruction [8], Linear Chirp Reconstruct  (LCR)
   [m_SR_Cl, m_LCR_Cl, IF_Cl] = R1_MR_and_LCR_grid(STFT, QM, Cs_simple, g, Lh, sigma_s, Nr, Nfft, N);

  
%         fprintf('VFB MB, ');
        [Cs_VFB_MB] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, 2, Nr);
%         Spl_MB = struct('spline', cell(1, Nr));
%         for m=1:Nr
%             Spl_MB(m).spline = spline((0:L-1)/L, (Cs_VFB_MB(m, :) - 1)*L/Nfft);
%         end
%         [m_SR_MB, m_LCR_MB, IF_MB] = R1_MR_and_LCR_spl(STFT, Spl_MB, g, Lh, sigma_s, Nr, Nfft, L);
        [m_SR_MB, m_LCR_MB, IF_MB] = R1_MR_and_LCR_grid(STFT, QM, Cs_VFB_MB, g, Lh, sigma_s, Nr, Nfft, N);




end

