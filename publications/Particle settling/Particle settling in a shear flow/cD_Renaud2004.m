function [c_D] = cD_Renaud2004( Re_p, n_PL )
%cD_Rp_Renaud2004 Power law drag law of Renaud et al. (2004)
%   Detailed explanation goes here

% Coefficient alpha
alpha_min = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);
alpha = 3.*(0.955)./(n_PL.^2+n_PL+1);
alpha_max = 3.*(0.955-0.021)./(n_PL.^2+n_PL+1);

% Coefficient X
X = alpha.*sqrt(6).^(n_PL+1)./6;

% Stokes drag
c_D0 = 24./Re_p.*X;

% Drag coefficient
c_D = c_D0 + 4.*0.44.*c_D0.^(11./24).*(0.558./(0.558+c_D0)).^(11./48)+0.44.*(0.558./(0.558+128.*c_D0)).^(11./12);

end % of function