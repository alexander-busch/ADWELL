function [ f ] = FanningFrictionFactor_Exp( rho, Usl, do, di, b, dpdl )
%FanningFrictionFactor_Exp Summary of this function goes here
%   Detailed explanation goes here

dh = HydraulicDiameter( do, di, b);

tauw = dpdl.*dh./4; % Wall shear stress
dynp = rho.*Usl.^2./2; % Dynamic pressure
f = tauw./dynp;

end

