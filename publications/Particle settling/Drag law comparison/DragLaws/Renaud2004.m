function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = Renaud2004( v_set, d_p, SR, mu_inf, mu_0, lambda_Cr, n_Cr, rho_f, rho_p )
%   Renaud2004 Renaud (2004)
%   Detailed explanation goes here



% Shear rate
SR_Standard = v_set./(d_p);

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR_Standard).^n_Cr);

% PL coefficients
[ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );

% Renaud (2004) Shear rate
alpha_min = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);
alpha = 3.*(0.955)./(n_PL.^2+n_PL+1);
alpha_max = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);
SR = alpha.*sqrt(6).*v_set./d_p;

% Particle Reynolds number 
Re_p  = rho_f.*v_set.^(2-n_PL).*d_p.^n_PL./K_PL;
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Renaud (2004) drag coefficient
X = alpha.*sqrt(6).^(n_PL+1)./6;
c_D0 = 24./Re_p.*X;
c_D = c_D0 + 4.*0.44.*c_D0.^(11./24).*(0.558./(0.558+c_D0)).^(11./48)+0.44.*(0.558./(0.558+128.*c_D0)).^(11./12);

% Settling velocity
v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;

end