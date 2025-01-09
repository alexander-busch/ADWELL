







SR = logspace(-2,4);



figure;
hold on;
set(gca,...
    'XScale','log',...
    'YScale','log');

k_Cr = 0.04119;
n_Cr =	0.4662;
mu_0 =	0.2058; 
mu_inf = 1.13e-08; 
eta = mu_inf+(mu_0-mu_inf)./((1+k_Cr.*SR).^n_Cr);
plot(SR,eta);


% Khatibi et al. (2016)
lambda_Cr = 0.0261;
n_Cr = 0.6082;
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);
plot(SR,eta);

% Busch et al. (2017)
lambda_Cr = 0.029;
n_Cr = 0.574;
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);
plot(SR,eta);




