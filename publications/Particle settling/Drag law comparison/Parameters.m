%% Fluid

% Density
rho_f = 1000;

% Rheological properties
load Cross.mat;
load Carreau.mat;

%% Solid

% Particle diameter [mm] --> [m]
d_p = [0.1; 1.16; 10] /1000;
% d_p = [1.16; ; 2.0 ; 3.0] /1000;
% d_p = [0.1; 0.3; 0.5; 1; 2; 3; 5; 10] /1000;

% Density
rho_p = 2560;

