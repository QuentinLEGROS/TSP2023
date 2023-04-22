function muc = compute_itSMAP(P,Ns,step_r,step_v)

[N,Nx] = size(P);
P1 = P;
muc = zeros(N,Ns);
% figure(11)
% imagesc(P1)
% pause(0.001)
for ns = 1:Ns
    [~,b] = max(P1(:)); % MAP
    mm=ceil(b/N); % time of max
%     nn = mod(b,N); % frequency of max
    nn = b-((mm-1)*N);
%     mm = (b/N-nn)*Nx;
    muc(nn,ns) = min(max(mm,1),Nx); % store MAP of current component
    P1(nn,min(max(mm-step_r,1),Nx):min(max(mm+step_r,1),Nx)) = 0; % discard associated vicinity in posterior
%     imagesc(P1)
%     pause(0.001)
    %forward
    for n = min(nn+1,N):N
%         disp('la')
        temp = min(max(mm-step_v,1),Nx):min(max(mm+step_v,1),Nx);
        [~,b] = max(P1(n,temp),[],2);
        muc(n,ns) = min(max(temp(b),1),Nx);
        P1(n,min(max(muc(n,ns)-step_r,1),Nx):min(max(muc(n,ns)+step_r,1),Nx)) = 0;
%         imagesc(P1)
%         pause(0.001)
        mm = temp(b);
    end
    %backward
    mm = muc(nn,ns);
    for n = max(nn-1,1):-1:1
        temp = min(max(mm-step_v,1),Nx):min(max(mm+step_v,1),Nx);
        [~,b] = max(P1(n,temp),[],2);
        muc(n,ns) = min(max(temp(b),1),Nx);
        P1(n,min(max(muc(n,ns)-step_r,1),Nx):min(max(muc(n,ns)+step_r,1),Nx)) = 0;
%         imagesc(P1)
%         pause(0.001)
        mm = temp(b);
    end
end
muc = min(max(muc,1),Nx);
[~,nn]=sort(mean(muc,1));
muc = muc(:,flip(nn));

% figure(5)
% subplot(2,1,1)
% imagesc(flipud(P'))
% subplot(2,1,2)
% plot(muc)
% ylim([1,Nx])
% pause(1)
