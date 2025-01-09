function [ f ] = FanningFrictionFactor_Ref( ReG, fluid, b  )
%FanningFrictionFactor_Ref Summary of this function goes here
%   Detailed explanation goes here

% Turbulent
t = ReG >= 2100 & strcmp(fluid, 'H2O');
f2 = t.*(0.3164./ReG.^0.25)./4; % Blasius (1912)

% Laminar
t = ~t;
f1 = t.*48./b./ReG;

f = f1 + f2;

end

