function [ val ] = g_func( m, M, L )
% [ val ] = g_func( m, M, L )
%
% compute the normalized and discretized square modulus of the Fourier transform of
% analysis window h
% 
% INPUT:
% s      : frequency bin(s) to compute
% M      : number of frequency bins to process
% L      : window time spread parameter
%
% OUTPUT:
% val    : g(m) values
%
% Author: D.Fourer (dominique.fourer@univ-evry.fr)
% Date: 23-feb-2021

C = (2*sqrt(pi)*L)/M;  %% normalizing constant
val = C * exp(-(2*pi*m/M).^2 * L^2);

end

