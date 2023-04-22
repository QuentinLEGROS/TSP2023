function Amp = Estim_Amp(muc,Y,w,G)

[N,~] = size(Y);
Ncomp = size(muc,2);
Amp = zeros(N,Ncomp);

for Nc = 1:Ncomp
    for n=1:N
        Amp(n,Nc) = (w(n,Nc)*sum(Y(n,:)))./G;
    end
end





