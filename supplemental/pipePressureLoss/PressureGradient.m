function [ dpdl ] = PressureGradient( f, rho, Usl, do, di, b)
%PressureGradient Summary of this function goes here
%   Detailed explanation goes here

dh = HydraulicDiameter( do, di, b);

dpdl = 2.*f.*rho.*Usl.^2./dh;
        
end

