%% this script is used to fit L with the function used by Nils for computing the window

clear
close all

addpath('tools');
addpath('synchrosqueezedSTFT');
addpath('Nils');

sigma_s = 0.09  %% the reference sigma_s for which L desired

N = 500;
M = 500;

%% Nils window
[g_nils,Lh] = create_gaussian_window(N, M, sigma_s);

er0 = inf;

%% find the best L matching sigma
for L_test = 2:1e-2:50
    %L = 18;
    L = L_test;
    A = 1/(sqrt(2*pi)*L);
    B = -1i * 2*pi / M;
    C = -1 / (2*L^2);
    gamma_K = 1e-4;
    K = 2 * L * sqrt(2*log(1/gamma_K));
    k = -K/2:K/2;
    k2 = k.^2;

    g_df = A * exp( C * k2);
    g_df = g_df/max(g_df);

    v1 = sum(g_nils(g_nils>1e-4));
    v2 = sum(g_df(g_df>1e-4));
    er = v2-v1;
    
    if abs(er) < abs(er0)
        L_best = L;
        er0 = er;
    end
    
end

L = L_best
    A = 1/(sqrt(2*pi)*L);
    B = -1i * 2*pi / M;
    C = -1 / (2*L^2);
    gamma_K = 1e-4;
    K = 2 * L * sqrt(2*log(1/gamma_K));
    k = -K/2:K/2;
    k2 = k.^2;

    g_df = A * exp( C * k2);
    g_df = g_df/max(g_df);

[~,I] = max(g_nils);
plot((-Lh):Lh, g_nils);

hold on
plot(k,g_df, 'r-.')


%% compare area under the curve:
sum(g_nils(g_nils>1e-4))
sum(g_df(g_df>1e-4))
