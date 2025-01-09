function [ dh ] = HydraulicDiameter( do, di, b)
%HydraulicDiameter Summary of this function goes here
%   Detailed explanation goes here

if b == 3 % Pipe
    dh = do;
else % Annulus
    dh = do - di;
end

end

