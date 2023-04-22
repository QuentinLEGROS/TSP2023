function [Mask_out,tempdata] = pseudoBay(tfr,Ncomp,M,L,div,beta,alpha,ds,Pnei,ifplot)
%
% Main algorithm
%
% INPUT:
% tfr           : Time-frequency represerntation of the MCS
% Ncomp         : Number of component
% M             : nombre de bins frequentiels
% L             : taille de la fenetre d'analyse en samples
% div           : Entropy choice (1 : KL | 2 : beta | 3 : Renyi)
% beta          : beta divergence hyperparameter
% alpha         : Renyi divergence hyperparameter
% ds            : variance of the random walk in the temporal model
% ifplot        : Boolean for debugging
% Pnei          : number of neighbors considered in the mask
%
%
% OUTPUT:
% Mask_out      : Mask of five bins per ridge. Modify code below to change
%                 the neighborhood



if ~exist('div', 'var')
 div=1;% default KL
end
if (~exist('beta', 'var') && (div==1))
 beta=1;% beta hyperparameter for beta divergence
end
if (~exist('alpha', 'var') && (div==2))
 alpha=0.5;%  alpha hyperparameter for Renyi divergence
end
if ~exist('ds', 'var')
 ds=3;% variance of the random walk in the temporal model
end
if ~exist('Pnei', 'var')
 Pnei=4;% default KL
end



data = transpose(abs(tfr(1:round(M/2),:))).^2; % Absolute value of the top half TFR

[Niter,N] = size(data); % Extract dimenssions



[F_mat, Fbeta_mat, IntFbeta, Falpha_mat, varF]=compFbeta_STFT(beta,alpha,M,L,N);
LFM = log(F_mat+eps);
% [tt,mmt]=max(data(100,:));
% [ttf,mmtf]=max(F_mat(100,:));
% Fm2 = sqrt(F_mat);
% 
% figure
% hold on
% plot(data(100,:),'b')
% plot(circshift(tt*(F_mat(100,:)./max(F_mat(100,:))),mmt-mmtf),'r')
% plot(circshift(tt*(Fm2(100,:)./max(Fm2(100,:))),mmt-mmtf),'k')
% legend('data','F','F2')


%% Initialization
tf=zeros(Niter,Ncomp); %Array to store the means of the depth
stf=zeros(Niter,Ncomp); %Array to store the variances of the depth 
tempdata = zeros(Niter,N,Ncomp);

for Nc = 1:Ncomp
    M2=floor(N/2); % mean of the depth prior when initizing (this can be changed)
    S2=N^2/12; % variance of the depth prior when initizing (this can be changed) % m and s2 close to the uniform distribution
 
    for t=1:Niter
        Y=data(t,:); % load current frame
        %% Main algorithm
        [M2,S2]=online_2D(Y+eps,LFM,Fbeta_mat,Falpha_mat,ds,M2,S2,IntFbeta,beta,alpha,div);
        %% Store values
        tf(t,Nc)=round(M2);
        stf(t,Nc)=S2;
    end
        
        %% Backward estimation
    for t=Niter:-1:1
        Y=data(t,:); % load current frame
        %% Main algorithm
        [M2,S2]=online_2D(Y+eps,LFM,Fbeta_mat,Falpha_mat,ds,M2,S2,IntFbeta,beta,alpha,div);
        %% Store values
        tf(t,Nc)=round(M2);
        stf(t,Nc)=S2;

    end

    % remove ridge using the three sigma rule of thumb
    % Computation of the ridge to remove
    for p = 1:Niter
        tempdata(p,:,Nc) = normpdf(1:N,tf(p,Nc),varF); % compute ridge vicinity
        tempdata(p,:,Nc) = data(p,tf(p,Nc)).*(tempdata(p,:,Nc)./max(tempdata(p,:,Nc))); % normalize
    end
    % Ridge removal
    
    
    
    if ifplot
        figure(2)
        subplot(Ncomp,2,(Nc-1)*2+1)     
        imagesc(transpose(data))
        yticklabels({'200','100','0'})
        yticks([50,150,250])
        title('Current data')
        subplot(Ncomp,2,(Nc-1)*2+2)
        imagesc(squeeze(transpose(tempdata(:,:,Nc))));title(strcat([num2str(Nc),'th estimated ridge']))
    end
    
    
    % Update the data without the just estimated rifge
    data = max(data - tempdata(:,:,Nc), 0);
end


% Computation of the mask
Mask_out=zeros(N*Niter,Ncomp);
veccol = transpose((1:Niter)-1).*N;
% for Nc = 1:Ncomp
%     Mask_out(max(tf(:,Nc)+veccol-1-5,ones(Niter,1)),Nc)=1;
%     Mask_out(max(tf(:,Nc)+veccol-1-4,ones(Niter,1)),Nc)=1;
%     Mask_out(max(tf(:,Nc)+veccol-1-3,ones(Niter,1)),Nc)=1;
%     Mask_out(max(tf(:,Nc)+veccol-1-2,ones(Niter,1)),Nc)=1;
%     Mask_out(max(tf(:,Nc)+veccol-1-1,ones(Niter,1)),Nc)=1;
%     Mask_out(max(tf(:,Nc)+veccol-1,ones(Niter,1)),Nc)=1;
%     Mask_out(min(tf(:,Nc)+veccol-1+1,(N*Niter)*ones(Niter,1)),Nc)=1;
%     Mask_out(min(tf(:,Nc)+veccol-1+2,(N*Niter)*ones(Niter,1)),Nc)=1; 
%     Mask_out(min(tf(:,Nc)+veccol-1+3,(N*Niter)*ones(Niter,1)),Nc)=1;
%     Mask_out(min(tf(:,Nc)+veccol-1+4,(N*Niter)*ones(Niter,1)),Nc)=1; 
%     Mask_out(min(tf(:,Nc)+veccol-1+5,(N*Niter)*ones(Niter,1)),Nc)=1; 
% end

for Nc = 1:Ncomp
    Mask_out(max(tf(:,Nc)+veccol-1,ones(Niter,1)),Nc)=1;
    for pn = 1:Pnei
        Mask_out(max(tf(:,Nc)+veccol-1-pn,ones(Niter,1)),Nc)=1;
        Mask_out(min(tf(:,Nc)+veccol-1+pn,(N*Niter)*ones(Niter,1)),Nc)=1; 
    end
end

Mask_out = reshape(Mask_out,[N,Niter,Ncomp]);
% size(tempdata)
% size(Mask_out)
% tempdata = transpose(sqrt(tempdata));


%% modif by DF
if Ncomp < 2
   Mask_out = sum(Mask_out,3);
   Mask_out = [Mask_out;Mask_out(end:-1:1,:)];
%    tempdata = [tempdata;tempdata(end:-1:1,:)];
else
   mask2 = zeros(2*N,Niter,Ncomp);
%    tempd2 = zeros(2*N,Niter,Ncomp);transpose(sqrt(tempdata));
   for i = 1:Ncomp
     mask2(:,:,i) = [Mask_out(:,:,i);Mask_out(end:-1:1,:,i)];
%      tempd2(:,:,i) = [tempdata(:,:,i);tempdata(end:-1:1,:,i)];
   end 
   Mask_out = mask2; clear mask2;
%    tempdata = tempd2; clear tempd2;
end



% tttf=abs(tfr).*Mask_out;
% 
% figure(4)
% subplot(1,3,1)
% imagesc((Mask_out))
% subplot(1,3,2)
% imagesc(abs(tfr))
% subplot(1,3,3)
% imagesc(tttf)
% 
% figure(5)
% subplot(1,3,1)
% plot(abs(tfr(:,200)))
% subplot(1,3,2)
% plot(tttf(:,200))
% subplot(1,3,3)
% plot(tempdata(:,200))
% pause(0.01)




