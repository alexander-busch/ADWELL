% Shah (2007) example

n_PL = 0.76;
K_PL = 6.51e-2;
rho_f = 1.0011e3;
d_p = 6.4e-4;
rho_p = 2.65e3;


% Shah (2007) coefficients
A = 6.9148.*n_PL.^2 - 24.838.*n_PL + 22.642;
B = -0.5067.*n_PL.^2 + 1.3234.*n_PL - 0.1744;

% Shah (2007) non-dimensional term
c_D_star = ((13.08.^(2-n_PL)./(2.^(2.*n_PL-2))) .* ...
    d_p.^(n_PL+2).*rho_f.^n_PL.*(rho_p-rho_f).^(2-n_PL)./K_PL.^2 ).^(0.5);

% Particle Reynolds number 
Re_p = (c_D_star./A).^(1./B);
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Settling velocity
v_set = (2.^(n_PL-1).*K_PL.*Re_p./(d_p.^n_PL.*rho_f)).^(1./(2-n_PL));

% Drag coefficient
c_D = 4.*d_p.*9.81./(3.*v_set.^2).*((rho_p./rho_f)-1);

















A = 6.9148.*n_PL.^2 - 24.838.*n_PL + 22.642;
B = -0.5067.*n_PL.^2 + 1.3234.*n_PL - 0.1744;

c_D_star = ((13.08.^(2-n_PL)./(2.^(2.*n_PL-2))) .* ...
    d_p.^(n_PL+2).*rho_f.^n_PL.*(rho_p-rho_f).^(2-n_PL)./K_PL.^2 ).^(0.5);

Re_p = (c_D_star./A).^(1./B);




%% Shah 1986 example

n_PL = 0.52;
K_PL = 0.01 ;
rho_f = 1.02;
d_p = 0.0275;
rho_p = 2.65;

% Shah (1986) coefficients
A = 21.394.*n_PL.^2 - 46.526.*n_PL + 29.899;
B = -2.9188.*n_PL.^4 + 6.7968.*n_PL.^3 - 5.6741.*n_PL.^2 + 2.4301.*n_PL + 0.0005;
C = 0.2323.*n_PL.^(-2.961);

% Shah (1986) non-dimensional term
c_D_star = ((3.5778.^(2-n_PL).*0.02615./(36.^(2.*n_PL-2))) .* ...
    d_p.^(n_PL+2).*rho_f.^n_PL.*(rho_p-rho_f).^(2-n_PL)./K_PL.^2 ).^(0.5);

% Particle Reynolds number 
Re_p = ((c_D_star-C)./A).^(1./B);
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Settling velocity
v_set = (36.^(n_PL-1).*K_PL.*Re_p./(0.1617.*d_p.^n_PL.*rho_f)).^(1./(2-n_PL));


% Conversion
n_PL = 0.52;
K_PL = 0.01 .* 4.448222e0 ./ 9.290304e-2;
rho_f = 1.02 * 1000;
d_p = 0.0275 / 39.37007874;
rho_p = 2.65 * 1000;













%% Debug

d_p = Cases.d_p(j,1);

mu_inf = Cross.mu_inf{k};
mu_0 = Cross.mu_0{k};
lambda_Cr = Cross.lambda{k};
n_Cr = Cross.n{k};
