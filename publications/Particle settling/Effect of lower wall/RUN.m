clear all;
close all;
clc;

%% Parameters
d_p = 0.002;
[Re_p,h]=meshgrid(logspace(-5,3),[0.01, 0.0075, 0.005, 0.0025, 0.001]);

% Ratio of particle radius to distance from the wall
ratio = d_p./(2.*h);

% Dimensionless wall distance
delta = (h-d_p./2)./d_p;

%% Create figure
fig = figure('Name','c_D = f(Re_p,h)','color','w');
hold on;
grid on;
box on;
title('Effect of wall distance on drag coefficient');
xlabel('Re_p [-]');
ylabel('h [mm]');
zlabel('c_D [-]');
set(gca,...
    'XScale','log',...
    'ZScale','log');

view([1,-1,1])

%% Standard drag coefficient (Schiller & Naumann 1933)
c_D0 = 24./Re_p;
c_D_SN = c_D0.*(1+0.15.*Re_p.^0.687);
% surf(Re_p,h*1000,c_D,'facecolor','none');
% plot(h,ones(1,length(h))*c_D);

%% Particle moving normal to the wall in quiescent fluid for Stokes flow, i.e. Re_p << 1 (Brenner 1961)
c_D = c_D0.*(1+ratio);
% surf(Re_p,h*1000,c_D./c_D_SN,'facecolor','none');
% plot(h,c_D);

%% Particle moving parallel to the wall in quiescent fluid (Faxen 1923), i.e. Re_p << 1
c_D = c_D0.*(1-9./16.*ratio+1./8.*ratio.^3-45./265.*ratio.^4-1./16.*ratio.^5).^(-1);
% surf(Re_p,h*1000,c_D./c_D_SN,'facecolor','none');
% plot(h,c_D);

%% Particle moving normal to the wall in shear flow (Zeng et al. 2009)
c_D0_s = 24./Re_p.*(1+0.138*exp(-2.*delta)+9./(16.*(1+2.*delta)));
alpha_s = 0.15-0.046.*(1-0.16.*delta.^2).*exp(-0.7.*delta);
beta_s = 0.687+0.066.*(1-0.76.*delta.^2).*exp(-delta.^0.9);
c_D_s = c_D0_s.*(1+alpha_s.*Re_p.^beta_s);
surf(Re_p,h*1000,c_D_s./c_D_SN,'EdgeColor','blue','facecolor','none');
% plot(h,c_D_s);

%% Particle moving parallel to the wall in quiescient fluid (Zeng et al. 2009)
c_D0_t = 24./Re_p.*(1.028-0.07./(1+4.*delta.^2)-8./15.*log(270.*delta./(135+256.*delta)));
alpha_t = 0.15.^(1-exp(-delta.^(0.5)));
beta_t = 0.687+0.313.*exp(-2.*delta.^(0.5));
c_D_t = c_D0_t.*(1+alpha_t.*Re_p.^beta_t);
surf(Re_p,h*1000,c_D_t./c_D_SN,'EdgeColor','green','facecolor','none');
% plot(h,c_D_t);


