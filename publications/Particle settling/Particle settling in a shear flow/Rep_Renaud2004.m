function [Re_p] = Rep_Renaud2004( rho_f, v_slip, n_PL, K_PL, d_p )
%   Rep_Renaud2004 Particle Reynolds number of Renaud et al. (2004)
%   Detailed explanation goes here

Re_p = rho_f.*v_slip.^(2-n_PL).*d_p.^n_PL./K_PL;

end % of function