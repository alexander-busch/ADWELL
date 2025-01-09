function [gamma_dot] = SR_Renaud2004( v_slip, n_PL, d_s )
%   Rep_Renaud2004 Particle Reynolds number of Renaud et al. (2004)
%   Detailed explanation goes here

% Coefficient alpha
alpha_min = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);
alpha = 3.*(0.955)./(n_PL.^2+n_PL+1);
alpha_max = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);

% Shear rate
gamma_dot = alpha.*sqrt(6).*v_slip./d_s;

end % of function