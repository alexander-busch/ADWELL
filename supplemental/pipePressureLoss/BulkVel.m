function [ bulkvel ] = BulkVel( mdot, rho, do, di, beta )
%BulkVel Summary of this function goes here
%   Detailed explanation goes here

if beta == 3 % Pipe
    di_help = 0;
else
    di_help = di;
end

VolFlowRate = mdot./rho;
XsectArea  = pi/4*(do^2-di_help^2);
bulkvel = VolFlowRate ./ XsectArea;

end

