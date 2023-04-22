function [gamma_e] = noise_level(STFT)
%NOISE_LEVEL Estimation of the noise level.
%   [gamma] = noise_level(STFT)
%
% INPUTS:
%   STFT    : Noisy short time fourier transform.
%
% OUTPUTS:
%   gamma_e : Estimation of the noise level, see [1].
%
% REFERENCES:
% [1] D. Donoho and I. Johnstone, “Ideal spatial adaptation via wavelet
% shrinkage,” Biometrika, vol. 81, pp. 425–455, 1994.

    gamma_e = median(abs(real(STFT(:))))/0.6745;
end