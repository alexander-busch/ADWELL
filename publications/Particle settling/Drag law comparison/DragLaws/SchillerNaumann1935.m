function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = SchillerNaumann1935( v_set, d_p, SR, mu_inf, mu_0, lambda_Cr, n_Cr, rho_f, rho_p )
%   SchillerNaumann1935 Schiller & Naumann (1935)
%   Detailed explanation goes here


% Shear rate
SR = v_set./(d_p);

% Apparent viscosity
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);

% Particle Reynolds number 
Re_p = Rep_Standard(rho_f, v_set, eta, d_p);
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Schiller & Naumann (1935) drag coefficient
c_D = (24./Re_p).*(1+0.15.*(Re_p.^0.687));

% Schiller & Naumann (1935) high Re correction
c_D(Re_p > 1000)=0.44;

% Settling velocity
v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;

end