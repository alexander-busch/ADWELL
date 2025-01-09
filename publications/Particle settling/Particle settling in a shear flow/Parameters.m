%% Domain

% Outer & inner diameter [inch] combinations acc. to K&M
d_o = [17.5; 12.25; 9.875; 8.5; 6.125];
d_i = [6.625; 6.625; 5; 4.5; 3.5];

% Hydraulic diameter [m]
d_h = (d_o - d_i) * 25.4/1000;

% X-Sectional are
A = pi./4.*((d_o.*25.4./1000).^2-(d_i.*25.4./1000).^2);

% Rotational speed of drill pipe
rpm = 150; % [50; 100; 150; 200];

%% Fluid

% Volumetric flow rate range [lpm] --> [m³/s]
VolFlowRate = [1; 10; 100; 1000; 10000] / 1000/60;

% Density
rho_f = 1000;

% Rheological properties
load CrossCoefficients.mat;

%% Solid

% Particle diameter [mm] --> [m]
d_p = [0.1; 1; 10] /1000;
% d_p = [0.1; 0.3; 0.5; 1; 2; 3; 5; 10] /1000;

% Density
rho_p = 2650;

