function [Re_p] = Rep_Standard( rho_f, v_slip, eta, d_p )
%   Rp_Standard Particle Reynolds number
%   Detailed explanation goes here

Re_p = rho_f.*v_slip.*d_p./eta;

end % of function

