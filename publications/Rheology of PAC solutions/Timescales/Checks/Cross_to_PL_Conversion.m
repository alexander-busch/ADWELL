fig1 = figure('color','w');
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'box','on');


gamma_dot = logspace(-2,4);

% Cross
eta = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*gamma_dot).^n_Cr);
plot(gamma_dot,eta);

gamma_dot_test = 100

% Metzner & Reed differential solution
eta = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*gamma_dot_test).^n_Cr);
tau = eta * gamma_dot_test;
epsilon = 1e-2*gamma_dot_test;
gamma_dot_up = gamma_dot_test+epsilon;
gamma_dot_down = gamma_dot_test-epsilon;
eta_up = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*gamma_dot_up).^n_Cr);
eta_down = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*gamma_dot_down).^n_Cr);
n_dash = (log(eta_up*gamma_dot_up)-log(eta_down*gamma_dot_down))/(log(gamma_dot_up)-log(gamma_dot_down))
K_dash = tau/gamma_dot_test^n_dash;
n_PL = n_dash;
k_PL = K_dash/((3.*n_dash+1)./(4.*n_dash))^n_dash;
plot(gamma_dot,k_PL.*gamma_dot.^(n_PL-1));


% Analytical solutuion
k_PL = exp(log(gamma_dot)*(ny_0-ny_inf)*n_Cr*k_Cr*gamma_dot/((ny_inf*(k_Cr*gamma_dot+1)^n_Cr+ny_0-ny_inf)*(k_Cr*gamma_dot+1)))*(ny_inf+(k_Cr*gamma_dot+1)^(-n_Cr)*ny_0-(k_Cr*gamma_dot+1)^(-n_Cr)*ny_inf)
n_PL=(-n_Cr*k_Cr*gamma_dot*ny_0+n_Cr*k_Cr*gamma_dot*ny_inf+k_Cr*gamma_dot*(k_Cr*gamma_dot+1)^n_Cr*ny_inf+k_Cr*gamma_dot*ny_0-k_Cr*gamma_dot*ny_inf+ny_inf*(k_Cr*gamma_dot+1)^n_Cr+ny_0-ny_inf)/(k_Cr*gamma_dot*(k_Cr*gamma_dot+1)^n_Cr*ny_inf+k_Cr*gamma_dot*ny_0-k_Cr*gamma_dot*ny_inf+ny_inf*(k_Cr*gamma_dot+1)^n_Cr+ny_0-ny_inf)
plot(gamma_dot,k_PL.*gamma_dot.^(n_PL-1));

legend('Cross fit (PAC4)',...
    'Metzner & Reed approximation',...
    'Analytical solution',...
    'location','northwest');

xlabel('Shear rate 1/s]');
ylabel('eta [Pa.s]');
