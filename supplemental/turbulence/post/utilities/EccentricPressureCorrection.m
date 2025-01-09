function [ dpdx_ec ] = EccentricPressureCorrection( dpdx_cc, Re, Re_cr, d_i, d_o, e, n  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R=zeros(length(Re),1);

for ii = 1:length(Re)
   if Re(ii)<Re_cr
        % Laminar
        R(ii) = 1 - 0.072*e/n*(d_i/d_o)^0.8454 - 1.5*e^2*n^0.5*(d_i/d_o)^0.1852 + 0.96*e^3*n^0.5*(d_i/d_o)^0.2527;
   else
        % Turbulent
        R(ii) = 1 - 0.048*e/n*(d_i/d_o)^0.8454 - 2/3*e^2*n^0.5*(d_i/d_o)^0.1852 + 0.285*e^3*n^0.5*(d_i/d_o)^0.2527;
   end
end
    
% % Laminar
% R = 1 - 0.072*e/n*(d_i/d_o)^0.8454 - 1.5*e^2*n^0.5*(d_i/d_o)^0.1852 + 0.96*e^3*n^0.5*(d_i/d_o)^0.2527;
% 
% % Turbulent
% R = 1 - 0.048*e/n*(d_i/d_o)^0.8454 - 2/3*e^2*n^0.5*(d_i/d_o)^0.1852 + 0.285*e^3*n^0.5*(d_i/d_o)^0.2527;

dpdx_ec = dpdx_cc.*R;


end

