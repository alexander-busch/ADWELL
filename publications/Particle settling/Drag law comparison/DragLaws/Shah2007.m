function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = Shah2007( v_set, d_p, SR, mu_inf, mu_0, lambda_Cr, n_Cr, rho_f, rho_p )
%   Shah2007 Shah (2007)
%   Detailed explanation goes here



% Shear rate
SR = v_set./(d_p);

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);

% PL coefficients
[ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );

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

end
