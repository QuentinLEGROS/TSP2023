function [X_hat, mask, Cs] = Brevdo_modeExtract_pp(tfr, L, Ncomp, K, step_r, step_v,bolpol)
%  [X_hat ] = Brevdo_modeExtract(tfr, Ncomp, lambda, clwin, K)
%
%  Mode retrieval and Mask extraction from an input time-frequency
%  representation
%
%  INPUT:
%  tfr : time-frequency representation
%  Ncomp : number of components to extract
%  K : vicinity of the ridge
%  lambda : ridge detection parameter 1 (regularization term) (default: 1e-1)
%  clwin :  ridge detection parameter 2
%
%  OUTPUT:
%
%
%
%  Author : Dominique fourer (dominique@fourer.fr)
%  Date : 12-feb-2021


%[tfr, stfr] = tfrsgab2(s, M, L, gamma_K);     %% (faster) classical  synchrosqueezed STFT


  

[M,N] = size(tfr);
Mh = round(M/2);
lambda = 0.0001; %1e-1;    %% ridge detection parameter 1
clwin  = 9;
%m_range(1:Mh)
[ Cs,mask ] = ridge_detect_brvmask_pp(abs(tfr(1:Mh,:)), (1:Mh)-1, Ncomp, lambda, clwin, K, step_r, step_v);


mask2 = zeros(M,N,Ncomp);
for i = 1:Ncomp
  mask2(:,:,i) = [mask(:,:,i);mask(end:-1:1,:,i)];
end
mask = mask2;


% %% compute the TF mask of each component (oracle computed from LM-reassigned spectrogram)
% mu        = 0.8;
% threshold = 0.2;
% %%[0.8 0.2] 17.84, 21.32, 22.21
% [ mask ] = compute_mask_comp(S, M, mu, threshold);


%% reconstruct the components
[ X_hat ] = sst_comp_ext(tfr, mask, L);


%% Postprocessing
if bolpol
    Cs = poly_int(Cs',step_r);
else
    Cs = Cs';
end



end

