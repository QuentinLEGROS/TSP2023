function [x] = istft(tfr,h,varargin)
%ISTFT inverse short time Fourier transform
%   [x] = istft(s,h)
%
% INPUTS:
%   tfr      : stft of x
%   h        : the filter
%
%   --  List of possible name-value pair argument
%   'cas'    : values 1, 2 and 3 respectively refers to
%              equations (3), (6) and (9) in paper [1]
%              default is 1.
%   'len'    : length of the signal to be reconstructed.
%              default is size(s, 2).
%   'shift'  : parameter used to shift the downsampling
%              default is 0.
%
% OUTPUTS:
%   x        : reconstruction from tfr
%
% REFERENCES:
% [1] S. Meignen and D.-H. Pham, “Retrieval of the modes of multi-
% component signals from downsampled short-time Fourier transform,”
% IEEE Trans. Signal Process., vol. 66, no. 23, pp. 6204–6215, Dec.
% 2018.
 
defaultCas = 1;
defaultLen = size(tfr, 2);
defaultShift = 0;

p = inputParser;
addRequired(p,'tfr');
addRequired(p,'h');
addParameter(p,'cas',defaultCas);
addParameter(p,'len',defaultLen);
addParameter(p,'shift',defaultShift);
parse(p,tfr,h,varargin{:});
cas = p.Results.cas;
Nsig = p.Results.len;
shift = p.Results.shift;

 [N,xrow_down] = size(tfr);
 downsamp = floor(Nsig/xrow_down);
 
 if (cas == 1)
  %case without periodizing
  x = zeros(Nsig,1);   
  Lh = (length(h)-1)/2;
  for icol=1:xrow_down,
   for q = 1:downsamp,
    ind = (icol-1)*downsamp+q; %mR+q
    if (ind <= Nsig)
     x(ind) = 1/h(Lh+q-shift)*mean(tfr(:,icol).*exp(2*pi*1i*(q-1-shift)*(0:N-1)'/N));
    end
   end
  end  
 end
 
 if (cas == 2)
  Lh = (length(h)-1)/2;
  x = zeros(Nsig,1);
  for i = 1:Nsig   
   ind  = 1+(floor(((i-1)-Lh-shift)/downsamp):ceil(((i-1)+Lh+shift)/downsamp));
   ind = ind(((i-1)-downsamp*(ind-1)-shift >= -Lh)&((i-1)-downsamp*(ind-1)-shift <= Lh));

   if (i > Lh)&&(i <= Nsig -Lh)
    x(i) = mean((tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           h(Lh+1+(i-1)-downsamp*(ind-1)-shift))/norm(h(Lh+1+(i-1)-downsamp*(ind-1)-shift))^2;
   else

    %To work well downsamp has to divide Nsig (cf paper)
    x(i) = mean((tfr(:,1+rem((ind-1)+xrow_down,xrow_down)).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))...
                *h(Lh+1+(i-1)-downsamp*(ind-1)-shift))/norm(h(Lh+1+(i-1)-downsamp*(ind-1)-shift))^2;     
   end 
  end
 end
 
 if (cas == 3)
  Lh = (length(h)-1)/2;
  x = zeros(Nsig,1);
  
  for i = 1:Nsig   
   ind  = 1+(floor(((i-1)-Lh-shift)/downsamp):ceil(((i-1)+Lh+shift)/downsamp));
   ind = ind(((i-1)-downsamp*(ind-1)-shift >= -Lh)&((i-1)-downsamp*(ind-1)-shift <= Lh));
   if (i > Lh)&&(i <= Nsig-Lh)
    x(i) = mean((tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           ones(length(Lh+1+(i-1)-downsamp*(ind-1)),1))/sum(h(Lh+1+(i-1)-downsamp*(ind-1)-shift));
   else
    x(i) = mean((tfr(:,1+rem((ind-1)+xrow_down,xrow_down)).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           ones(length(Lh+1+(i-1)-downsamp*(ind-1)),1))/sum(h(Lh+1+(i-1)-downsamp*(ind-1)-shift));     
   end 
  end 
 end
end
