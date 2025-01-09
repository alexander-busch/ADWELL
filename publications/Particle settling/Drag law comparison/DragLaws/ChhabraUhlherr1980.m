function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = ChhabraUhlherr1980( v_set, d_p, SR, mu_inf, mu_0, lambda_Ca, n_Ca, rho_f, rho_p )
%   ChhabraUhlherr1980 Chhabra & Uhlherr (1980
%   Detailed explanation goes here


% Shear rate
SR = v_set./(d_p);

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./((1+(lambda_Ca.*SR).^2).^((n_Ca-1)./2));

% Particle Reynolds number 
Re_p = Rep_Standard(rho_f, v_set, mu_0, d_p);
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Schiller & Naumann (1935) drag coefficient
c_D = (24./Re_p).*(1+0.15.*(Re_p.^0.687));

% Schiller & Naumann (1935) high Re correction
c_D(Re_p > 1000)=0.44;


% Chhabra & Uhlherr (1980) - Coefficients
A_1 = 0.65;
A_2 = 0.2;

% Chhabra & Uhlherr (1980) - Carreau number a.k.a. dimensionless time
Ca = 2.*lambda_Ca.*v_set./d_p;

% Chhabra & Uhlherr (1980) - X
X = 1+A_1.*(n_Ca-1).*Ca.^A_2;

% Chhabra & Uhlherr (1980) - "Final correlation of drag coefficient"
c_D = c_D.*X;

% Settling velocity
v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;

end

