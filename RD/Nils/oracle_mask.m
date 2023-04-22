function [mask] = oracle_mask(tfr, gamma)
% function [mask] = oracle_mask(tfr, gamma)
%
% Compute the mask corresponding to a given TFR
% considering values above a given threshold
%
% INPUT:
% tfr : input values
% gamma : threshold to use (default 1e-1)
%
%
% OUTPUT:
% mask: binary matrix of same dimension of tfr with values in {}
%
%
% Author: D.Fourer (dominique.fourer@univ-evry.fr)
% Date: 15-feb-2021

if ~exist('gamma', 'var')
  gamma = 1e-1;    
end

mask = (abs(tfr)>gamma);



