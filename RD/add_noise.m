function [y] = add_noise(x, n, SNR_in)
%ADD_NOISE add noise to x such that the input SNR is SNR_in
%   [y] = add_noise(x, n, SNR_in)
%
% INPUTS:
%   x        : signal
%   n        : noise
%   SNR_in   : input SNR
%
% OUTPUTS:
%   --  Output
%   y        : Noisy signal.

x = x(:);
n = n(:);

if SNR_in == inf
    y = x(:);
    return;
end

if SNR_in == -inf
    y = n(:);
    return;
end

Scale = sqrt(sum(abs(x).^2)/sum(abs(n).^2)*10^(-SNR_in/10));
y = x + Scale*n;
end

