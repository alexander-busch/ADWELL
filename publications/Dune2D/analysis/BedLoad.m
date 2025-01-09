
clear all;
clc;

fluid_bulk_vel = 0.3

gravity = 9.81

% Fluid properties
eta_f = 0.001002
rho_f = 1000
nu_f = eta_f/rho_f

% Solid properties
d_s = 0.0003
rho_s = 2650
bed_porosity = 0.45 

s = rho_s/rho_f

dstar = d_s*((s-1.0)*gravity/nu_f^2)^0.333333
teta_cr0 = (0.24/dstar)+0.055*(1-exp(-0.02*dstar))
teta_cr=teta_cr0;


stress_grav = rho_f * (s - 1) * gravity * d_s
teta = 0.416 / stress_grav

bed_porosity = 0.45
C_q_0 = 12/(1- bed_porosity)
q_0 = C_q_0 * sqrt(gravity * (s - 1.0) * d_s^3 * teta) * (teta - teta_cr)



% Shields number for inclined bed
phi = (23)/180*pi;

dy = -0.02:0.001:0.02;
dx = 0.02;

beta = atan(dy/dx);

theta_cr0 = 1;
figure, grid on, hold on;

% As used by Brors/Tron/...
theta_cr = theta_cr0.*sin(beta-phi)./sin(beta);
plot(beta./pi.*180,theta_cr);

% As given in Dey 
theta_cr = theta_cr0.*1./(cos(beta).*(1-tan(beta)./tan(phi)));
plot(beta./pi.*180,theta_cr);

% As used by Xiong(2014)
theta_cr = theta_cr0.*(cos(beta)-sin(beta)./tan(phi));
plot(beta./pi.*180,theta_cr);

% As used by Tron in bed_load.c - This is the correct formulation if beta
% is baed on the gradient as above.
theta_cr = theta_cr0.*sin(phi+beta)./sin(phi);
plot(beta./pi.*180,theta_cr);


