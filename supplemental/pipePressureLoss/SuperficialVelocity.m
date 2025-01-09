function [ Usl ] = SuperficialVelocity( mdot, rho, do, di, b )
%SuperficialVelocity Summary of this function goes here
%   Detailed explanation goes here

if b == 3 % Pipe
    di_help = 0;
else
    di_help = di;
end

VolFlowRate = mdot./rho;
XsectArea  = pi/4*(do^2-di_help^2);
Usl = VolFlowRate ./ XsectArea;

end

