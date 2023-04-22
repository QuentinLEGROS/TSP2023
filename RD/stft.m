function [tfr,norm2h] = stft(x,N,h,varargin)
%STFT short time Fourier transform
%   [TFR,norm2h] = stft(x,Nfft,h)
%
% INPUTS:
%   x        : signal
%   N        : number of frequency bins
%   h        : the filter
%
%   --  List of possible name-value pair argument
%   'cas'    : values 1, 2 and 3 respectively refers to
%              equations (3), (6) and (9) in paper [1]
%              default is 1.
%   'down'   : downsampling factor between 1 and floor(length(h)/2).
%              default is 1.
%   'shift'  : parameter used to shift the downsampling
%              default is 0.
% OUTPUTS:
%   --  Output
%   s        : the short time Fourier transform of x
%   norm2h   : the L2 norm of the filter on it support
%
% REFERENCES:
% [1] S. Meignen and D.-H. Pham, “Retrieval of the modes of multi-
% component signals from downsampled short-time Fourier transform,”
% IEEE Trans. Signal Process., vol. 66, no. 23, pp. 6204–6215, Dec.
% 2018.

defaultCas = 1;
defaultDown = 1;
defaultShift = 0;

p = inputParser;
addRequired(p,'x');
addRequired(p,'N');
addRequired(p,'h');
addParameter(p,'cas',defaultCas);
addParameter(p,'down',defaultDown); 
addParameter(p,'shift',defaultShift);
parse(p,x,N,h,varargin{:});
cas = p.Results.cas;
downsamp = p.Results.down;
shift = p.Results.shift;

Lh = floor(length(h)/2);
x = x(:);

if (length(h) >= N)
    error("window length too large : length(h) >= N");
end

if (downsamp > Lh)
    error("downsampling factor too large : 'down' > floor(length(h)/2)");
end

 
 [xrow,xcol] = size(x);
 
 t = 1+shift:downsamp:xrow; %the time instant, we consider the time instant shitfed by a factor shift.
  
 tfr= zeros (N,length(t)) ;
 if (cas == 1)
  %case without periodizing   
  trans  = zeros(1,length(t));
  norm2h = zeros(1,length(t));
  
  for icol=1:length(t),
   tau = -min([Lh,t(icol)-1]):min([Lh,xrow-t(icol)]); 
   tfr(1:length(tau),icol) = x(t(icol)+tau,1).*h(Lh+1+tau);
   trans(icol)  = tau(1);
   norm2h(icol) = norm(h(Lh+1+tau)); %we compute the L2 norm of the filter on its support 
  end
  tfr=fft(tfr,N); 
  A = exp(-2/N*pi*1i*(0:N-1)'*trans);
  tfr = tfr.*A;
 end
 
 if (cas == 2)||(cas == 3)
  %case with periodization    
  tau = -Lh:Lh;
  for icol = 1:length(t), 
   if (t(icol) > Lh) && (t(icol) <= xrow-Lh)
    tfr(1:length(tau),icol) = x(t(icol)+tau,1).*h(Lh+1+tau); 
   else
    tfr(1:length(tau),icol) = x(1+rem((t(icol)-1)+tau+xrow,xrow),1).*h(Lh+1+tau);
   end
   norm2h(icol) = norm(h(Lh+1+tau));%computation of the L2 norm of h 
  end
   tfr = fft(tfr,N);
   trans = Lh*ones(1,length(t)); 
   A = exp(2/N*pi*1i*(0:N-1)'*trans);
   tfr = tfr.*A;
 end
end 
