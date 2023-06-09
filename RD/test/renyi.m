function [RE_vec] = renyi(sigma_set, s_clean, L, Nfft)
    SL = length(sigma_set);
    RE_vec = zeros(1, SL);
    iSL = 0;
    alpha = 3;
    
    for sigma = sigma_set
        iSL = iSL + 1;
        fprintf('%u/%u\n', iSL, SL);
        [g, Lg] = gauss_win(L, sigma);
        [TFR, ~] = stft(s_clean, Nfft, g);
        Y = abs(TFR);
        TFR_MS = sum(Y(:));
        RE_vec(iSL) = 1/(1 - alpha)*log2(sum(Y(:).^alpha)/(TFR_MS^alpha)) - log2(L);
    end
    
end
