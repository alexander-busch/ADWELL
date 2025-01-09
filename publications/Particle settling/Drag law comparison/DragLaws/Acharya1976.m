function [ SR, eta, Re_p, Re_p_Standard, c_D, v_set ] = Acharya1976( v_set, d_p, SR, mu_Ca_inf, mu_Ca_0, lambda_Ca, n_Ca, rho_f, rho_p )
%   Acharya1976 Acharya et al. (1976)
%   Detailed explanation goes here

% Shear rate
SR = v_set./(d_p); % 

% PL coefficients based on current shear rate
if SR_p > (1.0/lambda_Ca) % Newtonian zero-shear viscosity region
    K_PL = mu_Ca_0;
	n_PL = 1.0;
else % Shear-thinning region
    K_PL = mu_Ca_0*lambda_Ca^(n_Ca-1.0);
	n_PL = n_Ca;
end

% Acharya et al. (1976) coefficients
f1 = 3.^((3.*n_PL-3)./2)...
    .* (33.*n_PL.^5-64.*n_PL.^4-11.*n_PL.^3+97.*n_PL.^2+16.*n_PL)...
    ./ (4.*n_PL.^2.*(n_PL+1).*(n_PL+2).*(2.*n_PL+1));
f2 = 10.5.*n_PL - 3.5;
f3 = 0.32.*n_PL + 0.13;

% Drag law coefficient correction of Kawase and Ulbrecht (1981) */
f1 = 3.0^((3.0*n_PL-3.0)/2.0)...
	* (-22.0*n_PL^2.0+29.0*n_PL+2.0)...
	/ (n_PL*(n_PL+2.0)*(2.0*n_PL+1.0));	
		
% Drag law coefficient correction of Kawase and Moo-Young (1986) */
f1 = 3.0^((3.0*n_PL-3.0)/2.0)...
	* (-7.0*n_PL^2.0+4.0*n_PL+26.0)...
	/ (5.0*(n_PL+2.0));

% Particle Reynolds number 
Re_p  = rho_f.*v_set.^(2-n_PL).*d_p.^n_PL./K_PL;
Re_p_Standard = Rep_Standard(rho_f, v_set, eta, d_p);

% Acharya et al. (1976) drag coefficient
if Re_p < 1
    c_D = 24./Re_p .* f1;   
else
    c_D = 24./Re_p .* f1 + f2./(Re_p.^f3);
end

% Settling velocity
v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;

end
