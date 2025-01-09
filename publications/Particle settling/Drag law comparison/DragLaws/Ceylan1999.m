function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = Ceylan1999( v_set, d_p, SR, mu_inf, mu_0, lambda_Cr, n_Cr, rho_f, rho_p )
%   Ceylan1999 Ceylan et al. (1999)
%   Detailed explanation goes here



% Shear rate
SR = v_set./(d_p); % 

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);

% PL coefficients
[ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );

% Particle Reynolds number 
Re_p  = rho_f.*v_set.^(2-n_PL).*d_p.^n_PL./K_PL;
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Ceylan et al. (1999) coefficients
c1 = 3.^(2.*n_PL-3)...
    .* (n_PL.^2-n_PL+3)...
    ./ (n_PL.^(3.*n_PL));
c2 = 4.*n_PL.^4./(24.*Re_p.^((n_PL-3)./3));

% Ceylan et al. (1999) drag coefficient
c_D = 24./Re_p .* (c1 + c2);   

% Settling velocity
v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;

end
