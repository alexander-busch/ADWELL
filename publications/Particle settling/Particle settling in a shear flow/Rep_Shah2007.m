function [Re_p] = Rep_Shah2007( rho_f, v_slip, n_PL, K_PL, d_s )
%   Rep_Shah2004 Particle Reynolds number of Shah et al. (2007)
%   Detailed explanation goes here

Re_p = rho_f.*v_slip.^(2-n_PL).*d_s.^n_PL./(2.^(n_PL-1).*K_PL);

end % of function