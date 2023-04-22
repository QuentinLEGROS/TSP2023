function tf = ChirprateRD(x,M,Ncomp,lambda_1,lambda_2,delta_1,delta_2)


N = length(x);
beta = - 400:400;
sigma = 0.2;

%% Gaussian window
[h, Lh] = create_gaussian_window(N,M,sigma);

%Time Frequency Chirprate Transform
[Cg_reassigned_store,~,~,~] = calcul_omega_beta(x,M,h, Lh,sigma,beta);


%% 3D RD algorithm
[tf] = Compute_IF(Cg_reassigned_store,lambda_1,lambda_2,delta_1,delta_2,Ncomp);


