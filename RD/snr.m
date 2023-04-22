function [snr_out] = snr(s,n)
%SNR signal to noise ratio
%   [snr_out] = SNR(s, n)
%
% INPUTS:
%   s        : signal
%   n        : noise
%
% OUTPUTS:
%   snr_out  : output SNR

if length(size(s)) == 2 && length(size(n)) == 2
    if min(size(s)) == 1 && min(size(n)) == 1
        s = s(:);
        n = n(:);
    end
end

rms = @(x) sqrt(mean(abs(x).^2));
snr_out = 20 * log10(rms(s)/rms(n));

end

