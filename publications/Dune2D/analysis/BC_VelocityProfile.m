
close all;
clear all;
clc;


H = 0.02; % channel_height
h = 0.0; % bed height
U_l = 0.085; % fluid bulk vel

A = pi()*H^2/4; % Pipe cross-sectional area
A_bed = (H/2)^2*acos(1-2*h/H)-(H/2-h)*(H*h-h^2)^(1/2); % Bed cross-sectional area
A_f = A-A_bed; % Flow x-sectional area of mixture 

Q = U_l*A; % Volumetric fluid flow rate

y = linspace(0,H,150);

figure;
hold on;

%% Laminar

% From pipe axis/channel half height to pipe radius/channel height
% Dimensional
u = 2.*U_l.*(1 - (y/H).^2); % h = radius, H = 0 (pipe axis)
plot(u,y);

% From pipe/channel lower wall to pipe/channel upper wall 
% Non-dimensional
y_nd = 2.*(y - 0.5*H) / H; 
u = 2.*U_l.*(1.0 - y_nd.^2);
plot(u,y);

% Dimensional
u=2.*U_l.*(1-(y-H./2).^2./(H./2).^2);  % h = channel y-coordinate from 0...H , H = channel height
plot(u,y);

U_2=Q/A_f;

u=2.*U_2-2.*U_2./((h-H)./2)^2.*(y-(h+(H-h)./2)).^2;
u=2.*U_2.*(1-1./((h-H)./2)^2.*(y-(h+(H-h)./2)).^2);

plot(u,y);

xlim([0 2*U_2]);
grid on;





%% Turbulent - Power law

m = 1/7; % exponent in turbulent flow eqns.

% PL - No bed height
uTmax = (m+1)*(m+2)*U_l/2; % max velocity for turbulent flow
plot(uTmax,h+(H-h)/2,'o');
u = uTmax*(2*(1-y/H)).^m; % turbulent velocity profile
plot(u,y);
u = uTmax*(2*(y/H)).^m; % turbulent velocity profile
plot(u,y);

% PL - With bed height
uTmax = (m+1)*(m+2)*U_2/2; % max velocity for turbulent flow
plot(uTmax,h+(H-h)/2,'o');
u = uTmax*(2*(1-(y-h)/(H-h))).^m; % turbulent velocity profile
plot(u,y);
u = uTmax*(2*((y-h)/(H-h))).^m; % turbulent velocity profile
plot(u,y);

% PL - With bed height

% UDF concept
u=zeros(1,length(y));
for i = 1:length(y)
    if y(i)<h
        u(i) = 0;
    elseif y(i)<(h+(H-h)/2)
        u(i) = uTmax*(2*((y(i)-h)/(H-h))).^m;
    else
        u(i) = uTmax*(2*(1-(y(i)-h)/(H-h))).^m;
    end
end
plot(u,y);


%% Turbulent - Log law
B = 2.5;

% Reynolds number
Re = rho_f * U_2 * d_h / eta_f;

% Friction factor (Haland formula, smooth pipe)
eps = 0;
lambda = ( 1/( -1.8*log(6.9/Re+(eps/3.7/d_h)^1.11) ) )^2;

% Wall shear stress
tau_w = 1/8*rho_f*lambda*U_2^2;

% Friction velocity
u_star = sqrt ( tau_w/rho_f );

% Mean velocity (nondimensional)
U_plus = U_2/u_star;

% Max velocity at pipe center
u_Tmax = (1 + 1.3258*sqrt(lambda))*U_2;

u_plus = u_Tmax + B*log(y/(d_h/2));











