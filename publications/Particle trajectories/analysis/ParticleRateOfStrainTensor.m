%% Clean up

clc;
clear all;
close all; 


%% Parameters

% Mean flow velocity
u = [0.0214 0.000188 0]';

% Shear rate of mean flow velocity
gamma_dot_m = 41.76;

% Particle diameter
d_p = 0.002;

% Particle velocity
v = [0.005 -0.007 0]';


% Particle slip velocity
v_slip = v-u;



%% Scalar formulation

% Particle slip velocity magnitude
v_slip_mag = sqrt((v_slip(1))^2 + (v_slip(2))^2 + (v_slip(3))^2);

% Mean flow velocity magnitude
u_mag = sqrt((u(1))^2 + (u(2))^2 + (u(3))^2);

% Shear rate based on particle slip velocity magnitude
gamma_dot_mag_sca = 2/d_p*v_slip_mag


%% Vector formulation

gamma_dot_vector = gamma_dot_m/u_mag*u+gamma_dot_mag_sca/v_slip_mag*v_slip;

gamma_dot_mag_vec = sqrt((gamma_dot_vector(1))^2 + (gamma_dot_vector(2))^2 + (gamma_dot_vector(3))^2);



%% Tensor formulation 

% Set diagonal to zero: dudx, dvdy, dwdz = 0
f = 1;

% Particle radius vector
r = d_p/2  * v_slip/v_slip_mag;

% Particle radius vector magnitude
r_mag = sqrt((r(1))^2 + (r(2))^2 + (r(3))^2);

dx = r(1);
dy = r(2);
dz = r(3);


% Particle deformation rate tensor
D(1,1) = f * v_slip(1) / dx;
D(2,1) = 1/2*(v_slip(1) / dx + v_slip(2) / dy);
D(3,1) = 1/2*(v_slip(1) / dz + v_slip(3) / dx);

D(1,2) = D(2,1);
D(2,2) = f * v_slip(2) / dy;
D(3,2) = 1/2*(v_slip(2) / dz + v_slip(3) / dy);

D(1,3) = D(3,1);
D(2,3) = D(3,2);
D(3,3) = f * v_slip(3) / dz;

% abs(sign(v_slip(i))
D

gamma_dot_mag_ten = sqrt((2*D(1,1))^2 + (2*D(2,2))^2 + (2*D(3,3))^2 + (D(2,1)+D(1,2))^2 + (D(3,1)+D(1,3))^2+ (D(2,3)+D(3,2))^2)

gamma_dot_mag_sca = 12/d_p*v_slip_mag