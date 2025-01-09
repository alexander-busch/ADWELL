function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = Shah1986( v_set, d_p, SR, mu_inf, mu_0, lambda_Cr, n_Cr, rho_f, rho_p )
%   Shah1986 Shah (1986)
%   Detailed explanation goes here

% Conversions






% Shear rate
SR = v_set./(0.5.*d_p);
% SR = 36.*v_set./d_p;

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);

% PL coefficients
[ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );
K_PL = K_PL ./ 4.448222e0 .* 9.290304e-2; % N/m² --> lbf / ft²

% Shah (1986) coefficients
A = 21.394.*n_PL.^2 - 46.526.*n_PL + 29.899;
B = -2.9188.*n_PL.^4 + 6.7968.*n_PL.^3 - 5.6741.*n_PL.^2 + 2.4301.*n_PL + 0.0005;
C = 0.2323.*n_PL.^(-2.961);

% Shah (1986) non-dimensional term, d_p converted to inch, K_PL converted to specific
% gravities
d_p = d_p * 39.37007874; %  m --> in
rho_f = rho_f/1000;
rho_p = rho_p/1000;
c_D_star = ((3.5778.^(2-n_PL).*0.02615./(36.^(2.*n_PL-2))) .* ...
    d_p.^(n_PL+2).*rho_f.^n_PL.*(rho_p-rho_f).^(2-n_PL)./K_PL.^2 ).^(0.5);

% Particle Reynolds number 
Re_p = ((c_D_star-C)./A).^(1./B);
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Settling velocity
v_set = (36.^(n_PL-1).*K_PL.*Re_p./(0.1617.*d_p.^n_PL.*rho_f)).^(1./(2-n_PL));

v_set =  v_set ./ 3.280839895; %  ft --> m

c_D = c_D_star;

end
