function result=Gamma(N,Fs,Delta,mu,psi,w,t)
  L=size(Delta,2);
  Mw=repmat(w,L,1);
%   t=time_axis(N,Fs);
  result=transpose(sum(Mw.*exp(mu.'*t+i*(Delta.'*t+psi.'*(t.^2)/2)),2));
  
  if(length(find(result==0))==0)
     % disp 'zeros should never happen...'
     return;
  end;
  %assert(length(find(result==0))==0);  % zeros should never happen...
  idx=[find(isnan(result)) find(isinf(result))];  % in case of problems with large mu value...
  result(idx)=1;  % ...reset to default: no correction

