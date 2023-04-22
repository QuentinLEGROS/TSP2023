function [tfr,Fc] = tfrCRgab2(x, M, L, gamma_K,beta)
% [tfr] = tfrgab(x, M, L, gamma_K)
% Compute the discrete-time chirprate Transform (using FFT)
% 
% INPUT:
% x      : the signal to process
% M      : number of frequency bins to process (default: length(x))
% L      : window duration parameter:  w0 * T, (default: 10)
% gamma_K: threshold applied on window values (default: 10^(-4))
%
% OUTPUT:
% tfr    : discrete Stockwellogram
%
% Author: Q.Legros
% Date: 22-01-2023
% Ref: [ZHU, Xiangxiang, YANG, Haizhao, ZHANG, Zhuosheng, /et al./Frequency-chirprate reassignment. /Digital Signal Processing/, 2020,vol. 104, p. 102783]

x = x(:).';          %convert x as a row vector
N = length(x);

if ~exist('M', 'var')
 M = N;
end
if ~exist('L', 'var')
 L = 10;
end
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end

tfr = zeros(M, N);
step = 1;
nb_beta = 10;
Fc = zeros(M, N, nb_beta);


K = 2 * L * sqrt(2*log(1/gamma_K));  %% window length in samples

A = 1/(sqrt(2*pi)*L);
B = -1i * 2*pi / M;
C = -1 / (2*L^2);
D = -1i * 2*pi / (2*M);

mm = m_axis(M);


for n = 1:N
  count = 1;
  k_min = min(n-1, round(K/2));
  k_max = min(N-n, round(K/2));
  k = (-k_min):k_max;
  k2 = k.^2;
  g = A * exp( C * k2);
  
  for bet = 0:step:(step*nb_beta)-step
    Fc(:,n,count) = fft(x(n+k) .* g, M) .* exp(B * mm * (n-1 - k_min)) .* exp(D  * bet * ((n-1 - k_min).^2));
    count = count+1;
  end
end
tfr = Fc(:,:,1);
end
