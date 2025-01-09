close all;
clear all;
clc;

SR = logspace(-5,4); % Shear rate range

C_P = 2560*9.81*0.01; % Solid pressure 


% Coefficients 
angle = 30.0; % Angle of repose [°] 
cohesion = 1.0; % Cohesion coefficient [Pa] 
tau0 = cohesion + C_P*tan(angle*pi/180.0); % Yield stress [Pa] 
k = 2.0; % Consistency index 
n = 0.5; % Flow index  
eta0 = 1.0e6; % Limiting dynamic viscosity for no flow 

% Strain rates
S_dot = 40;
S_crit = max(tau0, 1.0E-10) / max(eta0, 1.0E-10);


figure; hold on; set(gca,'XScale','log','YScale','log');

for i=1:length(SR)
    if (SR(i) <= S_crit)
        plot(SR(i),eta0,'k*');
    else
        % plot(SR(i),(tau0 + k*SR(i)^n - S_crit^n)/SR(i),'r*');
        plot(SR(i),(tau0 + k*SR(i)^n)/SR(i),'b*');
    end
end