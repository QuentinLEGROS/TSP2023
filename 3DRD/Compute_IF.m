function [tf] = Compute_IF(Fc,lambda_1,lambda_2,delta_r,delta_c,K)
% 
% Compute Algo 1 of the paper 'Frequency-chirprate reassignment'
% 
% %
% INPUTS:
% Fc: TFCR representation
% lambda_1: Constraint weight on omega
% lambda_2: Constraint weight on beta
% delta_r: allowable variations for r
% delta_c: allowable variations for c
% K : number of components
%
% OUTPUTS:
% tf  : estimated instentaneous frequency


[N,M,L]=size(Fc);
tf = zeros(N,K);
r = zeros(N,1);
c = zeros(N,1);

for k=1:K % loop over number of components
       % random initialization
       q = N/2;
       Fcc = squeeze(Fc(q,:,:));
       maximum = max(max(Fcc));
       ma=find(Fcc==maximum);
       % verify if there are not identical maxima
       if length(ma)>1
           rr = randperm(2);
           c(q)=floor(ma(rr(1))/M)+1;
           r(q) = ma(rr(1))-(M*floor(ma(rr(1))/M));
       else
           [r(q),c(q)]=find(Fcc==maximum);
       end
       
       % Forward
       for n = q+1:N
           % vicinity
           I_r = max(1,r(n-1)-delta_r):min(r(n-1)+delta_r,M);
           I_c = max(1,c(n-1)-delta_c):min(c(n-1)+delta_c,L);
           % compute maximum
           temp = squeeze(Fc(n,I_r,I_c)) - lambda_1*(I_r'-r(n-1)).^2 - lambda_2*(I_c-c(n-1)).^2;
           maximum = max(max(temp));
           % finding next value
           [r1,c1]=find(temp==maximum);
           r(n) = I_r(r1);
           c(n) = I_c(c1);
       end

       % backward
       for n = q-1:-1:1
           % vicinity
           I_r = max(1,r(n+1)-delta_r):min(r(n+1)+delta_r,M);
           I_c = max(1,c(n+1)-delta_c):min(c(n+1)+delta_c,L);
           % compute maximum
           temp = squeeze(Fc(n,I_r,I_c)) - lambda_1*(I_r'-r(n+1)).^2 - lambda_2*(I_c-c(n+1)).^2;
           maximum = max(max(temp));
           % finding previous value
           [r1,c1]=find(temp==maximum);
           r(n) = I_r(r1);
           c(n) = I_c(c1);
       end
       
       % store instentaneous frequency
       tf(:,k)=r;

        % Update the CT transform (to avoid estimating the same mode several times)
        for n=1:N
           Fc(n,max(1,r(n)-delta_r):min(r(n)+delta_r,M),max(1,c(n)-delta_c):min(c(n)+delta_c,L)) = 0;
        end
end

