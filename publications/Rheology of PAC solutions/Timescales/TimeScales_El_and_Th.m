%% Clean-up
clear all;
close all;
clc;

% Parameters

rho_f = 1000;

eta = logspace(-2,1);
l = logspace(-4,4);
[ETA,L] = meshgrid(eta,l);



% Rheological time scales
load Timescales_3ITT_100_50.mat;
load Timescales_3ITT_300_600_Large.mat;
load Timescales_FS.mat;
lambda_el_lower = min(MaxwellRelaxationTime);
lambda_el_upper = max(MaxwellRelaxationTime);
lambda_th_lower = (lambda_3ITT_100_50_Lab1+lambda_3ITT_100_50_Lab2)/2;
lambda_th_upper = (lambda_3ITT_300_600_Lab1_HighShearRate+lambda_3ITT_300_600_Lab2_HighShearRate)/2;



% Particle scale, Elasticity

lambda_el_upper = 1e-1;
El_upper = ETA.*lambda_el_upper./rho_f./L.^2;
lambda_el_lower = 1e-4;
El_lower = ETA.*lambda_el_lower./rho_f./L.^2;

fig_name = 'Particle scale, Elasticity';

fig_color = 'w';
fig_units = 'centimeters';
fig_position = [1 1 21 14.8]; % DIN A5 

fig_xlabel = 'L [m]';
fig_ylabel = '\eta [Pa.s]';
fig_zlabel = 'El [-]';

fig_xlim = [min(l) max(l)];
fig_ylim = [min(eta) max(eta)];
fig_zlim = [0.1 max(max(max(El_lower)),max(max(El_upper)))];

fig_scale = 'log';

Create_Figure3D(fig_name, fig_color, fig_units, fig_position, fig_xlabel, fig_ylabel, fig_zlabel, fig_xlim, fig_ylim, fig_zlim, fig_scale);
surf(L,ETA,El_upper,'FaceColor','[1 .5 0]');
surf(L,ETA,El_lower,'FaceColor','red');



% Particle scale, Thixotropy

lambda_th_upper = 5e3;
Th_upper = ETA.*lambda_th_upper./rho_f./L.^2;
lambda_th_lower = 1e1;
Th_lower = ETA.*lambda_th_lower./rho_f./L.^2;

fig_name = 'Particle scale, Thixotropy';
fig_zlabel = 'Th [-]';
fig_zlim = [0.1 max(max(max(El_lower)),max(max(El_upper)))];

Create_Figure3D(fig_name, fig_color, fig_units, fig_position, fig_xlabel, fig_ylabel, fig_zlabel, fig_xlim, fig_ylim, fig_zlim, fig_scale);
surf(L,ETA,Th_upper,'FaceColor','[1 .5 0]');
surf(L,ETA,Th_lower,'FaceColor','red');
