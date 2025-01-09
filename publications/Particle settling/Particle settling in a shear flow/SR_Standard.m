function [gamma_dot] = SR_Standard( v_slip, d_s )
%   Rep_Renaud2004 Particle Reynolds number of Renaud et al. (2004)
%   Detailed explanation goes here

gamma_dot = 2.*v_slip./d_s;

end % of function