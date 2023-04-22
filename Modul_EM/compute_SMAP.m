function muc = compute_SMAP(P,Ns,Mu)

[N,Nx] = size(P);
P1 = P;

for ns = 1:Ns
    [~,muc(:,ns)] = max(P1,[],2); % MAP
    for n = 1:N
%         P1(n,muc(n,ns)-10:muc(n,ns)+10) = -1;
        P1(n,min(max(muc(n,ns)-30,1),Nx):min(max(muc(n,ns)+30,1),Nx)) = 0;

    end
end


