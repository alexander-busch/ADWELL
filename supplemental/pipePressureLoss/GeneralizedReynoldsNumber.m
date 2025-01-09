function [ ReG ] = GeneralizedReynoldsNumber( rho, Usl, do, di, b, n, K )
%GeneralizedReynoldsNumber Summary of this function goes here
%   Detailed explanation goes here

dh = HydraulicDiameter( do, di, b);

Nom = rho.*Usl.^(2-n).*dh.^n;
DeN = (24./b).^(n-1) .* K .* ((1+n.*b)./(n+n.*b)).^n;
ReG = Nom ./ DeN;

end

