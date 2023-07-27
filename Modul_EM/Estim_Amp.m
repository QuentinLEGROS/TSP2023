function Amp = Estim_Amp(muc,Y,w,G,M,L)

[N,~] = size(Y);
Ncomp = size(muc,2);
Amp = zeros(N,Ncomp);

for Nc = 1:Ncomp
    for n=1:N
        Amp(n,Nc) = (w(n,Nc)*sum(Y(n,:)))./G;
    end
end

% sqrt(Amp(100,1))/2

% G2 = M/(sqrt(pi)*L);
% for Nc = 1:Ncomp
%     for n=1:N
%         Amp(n,Nc) = (w(n,Nc)*sum(Y(n,:)))./(G2);
%     end
% end

% n=100;
% G2 = M/(2*sqrt(pi)*L);
% Amp2 = (w(n,1)*sum(Y(n,:)))./(G2);
% Amp2 = sqrt(Amp2)
% 
% G3 = M/(2*sqrt(2)*L);
% Amp3 = (w(n,1)*sum(Y(n,:)))./G3;
% Amp3 = sqrt(Amp3)
% pause