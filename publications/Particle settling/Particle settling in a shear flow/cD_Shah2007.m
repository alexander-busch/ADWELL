function [c_D] = cD_Shah2007( Re_p, n_PL )
%cD_Shah2007 Power law drag law of Shah et al. (2007)
%   Detailed explanation goes here

% Coefficients A & B
A = 6.9148.*n_PL.^2 - 24.838.*n_PL + 22.642;
B = -0.5067.*n_PL.^2 + 1.3234.*n_PL - 0.1744;

% Drag coefficient
c_D = (A.^2.*Re_p.^(B)).^(1./(2-n_PL));

end % of function