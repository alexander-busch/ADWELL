function [ TauW ] = WallShearStress( do, dpdl )
%WallShearStress Summary of this function goes here
%   Detailed explanation goes here

TauW = do.*dpdl./4;

end

