function [Cg_reassigned_store,Cg_store,beta_store,omega_store] = calcul_omega_beta(s, Nfft, g, Lg, sigma_s,beta)
%
% INPUTS:
% s: real or complex signal
%    has to be a single mode
%    the phase of the signal has to be 0 at time t=0
% Nfft: frequency bins
% g: gaussian window s.t. g(x) = exp(-pi*x^2/sigma_s^2)
% Lg: gaussian window length
% sigma_s: gaussian coefficient
% beta : the set of potential chirp rates 
%
% OUTPUTS:
% beta_store  : the tensor of the values for beta
% omega_store : the tensor of the values for omega

s = s(:);

L = length(s);

t0  = ((0:2*Lg)'-Lg)/L;
g_chirp = g.*exp(-pi*1i*(t0.^2)*beta);
len_b = length(beta);

beta_store     = zeros(L,Nfft,len_b);
omega_store    = zeros(L,Nfft,len_b); 
Cg_store       = zeros(L,Nfft,len_b);
Cg_reassigned_store  = zeros(L,Nfft,len_b);

gamma = 10^(-6);

for n=1:L
%  n
 % Cg, window g
 time_inst = -min([Lg,n-1]):min([Lg,L-n]);% truncating Gaussian to remain on the signal support
 Cg= fft(s(n+time_inst).*g_chirp(Lg+time_inst+1,:),Nfft)/L;
    
 % C_xg, window xg
 C_xg = fft(s(n+time_inst).*((time_inst)'/L).*g_chirp(Lg+time_inst+1,:),Nfft)/L;
    
 % C_xxg, window xxg
 C_xxg = fft(s(n+time_inst).*(((time_inst)'/L).^2).*g_chirp(Lg+time_inst+1,:),Nfft)/L;

 % operator beta hat 
 denom_beta = C_xg.^2 - C_xxg.* Cg;
 beta_hat   = 1/(2*pi)*imag(Cg.^2./denom_beta) + ones(Nfft,1)*beta;

 
 % operator omega hat
 % we need a rounded value of beta hat to compute omega hat 
 beta_hat_round = 1+max(0,min(len_b-1,round(beta_hat-beta(1)))); %indices corresponding to beta_hat

 C_g_opt  = zeros(Nfft,len_b);
 C_xg_opt = zeros(Nfft,len_b); 
    
 for k=1:len_b
  for p = 1:Nfft
   C_g_opt(p,k)  = Cg(p,beta_hat_round(p,k));
   C_xg_opt(p,k) = C_xg(p,beta_hat_round(p,k));
  end
 end
 
 omega_hat = 1/(2*pi*sigma_s^2)*imag(C_xg_opt./C_g_opt)+ ((0:Nfft-1)*L/Nfft)';
 
 T = abs(Cg).^2;
 Cg_store(n,:,:) = T;
  
 % Storing beta hat
 beta_store(n,:,:) = beta_hat;
    
 %Storing omega hat 
 omega_store(n,:,:) = omega_hat;
 
 %% Reassignment step
 Cg_reassigned = zeros(Nfft,len_b);
 A = denom_beta > gamma^2; 
 B = T > gamma^2;
 C = round(1+omega_hat);
 
 for k = 1:len_b
  for p = 1:Nfft
   if (B(p,k) > gamma^2) && (A(p,k) > gamma^2)
    p1 = C(p,k);
    k1 = beta_hat_round(p,k);
    if (p1 >=1) && (p1 <= Nfft)
     Cg_reassigned(p1,k1)  = Cg_reassigned(p1,k1) + T(p,k);
    end
   end
  end
 end
 Cg_reassigned_store(n,:,:) = Cg_reassigned; 
end
 
 