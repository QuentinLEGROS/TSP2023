function wc = Mod_NR_multi(y,wc,Ft,alpha,Ns,muc,p)



% Perform second order gradient descent
% 
% INPUT:
% y           : Spectrogram at time n
% wc          : current mixture weight at time n
% p           : posterior distribution at time n
% Ft          : data distribution
% alpha       : Dirichlet prior parameters

% 
% OUTPUT:
% wc         : estimated mixtre weight
%
% Author: Q.Legros


%% parameters and matrix pre computation
err=10;
it=0;
[M,~,~] = size(Ft);
wc = transpose(wc);
H = zeros(length(wc),length(wc));
graf = zeros(length(wc),1);
%% Main loop

while (err>0.0001 && it<30 )
    w10=wc;       
    A = (1-sum(wc))/M;
    for ns = 1:Ns
        A = A + wc(ns)*Ft(:,muc(ns));
    end
    for ns = 1:Ns
       P(:,ns) = (Ft(:,muc(ns))-(1/M)) ./ A;
       graf(ns) = sum(y.*P(:,ns));
    end

    for ns = 1:Ns
        for ns2 = 1:Ns
            H(ns,ns2) = -sum(y.*(P(:,ns).*P(:,ns2)));
        end
    end

    % NR step
    wc = wc - pinv(H)*graf;
    wc = min(max(wc,1e-3),1-1e-3);
    if sum(wc) > 1
        wc = wc./((1+1e-3)*sum(wc));
    end
    %% Compute error
    err = sum(abs(wc-w10)./min(wc,w10));
    it = it + 1;

end

