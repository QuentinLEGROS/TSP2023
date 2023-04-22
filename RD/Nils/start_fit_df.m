%% this script is used to fit L with the function used by Nils for computing the window

clear
close all

addpath('tools');
addpath('synchrosqueezedSTFT');
addpath('Nils');

L = 20   %%reference L value
%sigma_s = 0.09
N = 500;
M = 500;

%% Df window
A = 1/(sqrt(2*pi)*L);
B = -1i * 2*pi / M;
C = -1 / (2*L^2);
gamma_K = 1e-4;
K = 2 * L * sqrt(2*log(1/gamma_K));
k = -K/2:K/2;
k2 = k.^2;

g_df = A * exp( C * k2);
g_df = g_df/max(g_df);
    
    

er0 = inf;

%% find the best matching sigma
for S_test = 0:1e-3:sigma_s
    sigma_s = S_test;
    [g_nils,Lh] = create_gaussian_window(N, M, sigma_s);

    v1 = sum(g_nils(g_nils>1e-4));
    v2 = sum(g_df(g_df>1e-4));
    er = v2-v1;
    
    if abs(er) < abs(er0)
        S_best = sigma_s;
    end
    er0 = er;
end

sigma_s = S_best
[g_nils,Lh] = create_gaussian_window(N, M, sigma_s);

[~,I] = max(g_nils);
plot((-Lh):Lh, g_nils);

hold on
plot(k,g_df, 'r-.')

%% compare area under the curve:
sum(g_nils(g_nils>1e-4))
sum(g_df(g_df>1e-4))
