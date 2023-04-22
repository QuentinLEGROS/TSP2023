function [a]=Compute_A_DF(x,Fs,m,L,n)
% function [a,mu,phi,omega,psi,delta_t]=reassignment(x,Fs,m)
% Estimate sine parameters at position m from fft(x)
%
% INPUT :
% x : signal
% Fs : sampling rate
% m : frequency index
%
% OUTPUT :
% a : amplitude at t=0 

N=length(x);
%% Compute time
%   gamma_K = 10^(-4);
%   K = 2 * L * sqrt(2*log(1/gamma_K));  %% window length in samples
%   A = 1/(sqrt(2*pi)*L);
%   B = -1 / (L^2);
%   C = -1 / (2*L^2);
% 
%   k_min = min(n-1, round(K/2));
%   k_max = min(N-n, round(K/2));
%   k = (-k_min):k_max;
%   k2 = k.^2;
%   w = A * exp( C * k2);
%   wd = B*k .* w;
%   wdd = ((B*k).^2) .* w;
% 
%   
%   x = x(k+n);
%   t = k;
%% Compute Amplitude for the IF associated to time

  % prepare the windows (using the Hann window, two-time differentiable)
  t=time_axis(N,Fs);

  w=gabor_window( N, L, n );
  wd=gabor_window_d1( N, L, n );
  wdd=gabor_window_d2( N, L, n );
  
  tw=t.*w;
  twd=t.*wd;
  % compute the spectra
  Xw=fft(zerophase_signal(w.*x));
  Xwd=fft(zerophase_signal(wd.*x));
  Xtw=fft(zerophase_signal(tw.*x));
  Xwdd=fft(zerophase_signal(wdd.*x));
  Xtwd=fft(zerophase_signal(twd.*x));
  % search for peaks
  if nargin<3
    % consider the global maximum
    [~,m]=max(abs(Xw));
  end
  % frequency
  delta_omega=-imag(Xwd(m)./Xw(m));
  % amplitude modulation
  mu=-real(Xwd(m)./Xw(m));
  % time offset
  % frequency modulation
  psi=(imag(Xwdd(m)./Xw(m))-imag((Xwd(m)./Xw(m)).^2))./(real((Xtw(m).*Xwd(m))/(Xw(m)).^2)-real(Xtwd(m)./Xw(m)));
  % then amplitude and phase
  p=Gamma(N,Fs,delta_omega,mu,psi,w,t);
  a=abs(Xw(m)./p);
  
  %%%% A'(omega)
end
