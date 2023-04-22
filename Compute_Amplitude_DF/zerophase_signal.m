function dst=zerophase_signal(src)
  N=length(src);
  if isodd(N)
    H=(N-1)/2;
  else
    H=N/2;
  end
  dst=shift(src,-H);
  % and +H to invert the zerophasing (+H=-H mod N when N is even)
end

function y=shift(x,b,dim)
  if (nargin<3)||(dim==2)
    y=(circshift(x.',b)).';
  else
    y=(circshift(x,b));
  end
end
