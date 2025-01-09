%% Clean-up
clear all;
close all;
clc;

addpath 'E:\OneDrive - NTNU\SWwd\MATLAB\Generic\m-files';



% Load time scales
load Timescales_3ITT_100_50.mat;
load Timescales_3ITT_300_600.mat;
load Timescales_3ITT_300_600_Large.mat;
load Timescales_FC.mat;
load Timescales_FNSC.mat;
load Timescales_FS.mat;


% Parameters
global rho_f rho_s k_Cr n_Cr ny_0 ny_inf;
global lambda_Maxwell_1 lambda_Maxwell_2 lambda_3ITT_1 lambda_3ITT_2;
global fig_MFS_3ITT fig_MFS_Maxwell fig_PS_3ITT fig_PS_Maxwell colorlist;

% Outer & inner diameter [inch] combinations acc. to K&M
d_o = [17.5; 12.25; 9.875; 8.5; 6.125];
d_i = [6.625; 6.625; 5; 4.5; 3.5];

% Hydraulic diameter [meter]
d_h = (d_o - d_i) * 25.4/1000;

% Length scale of Main Flow Scale
L = 1; % [1; 10; 100; 1000];

% Volumetric flow rate range [lpm] --> [m³/s]
VolFlowRate = [1; 10; 100; 10000] / 1000/60;

% Max. volumetric flow rate as f(d_h) [lpm] --> [m³/s] acc. to K&M
VolFlowRateMax = [5700; 4200; 3500; 2300; 700] / 1000/60;


% Densities
rho_f = 1000;
rho_s = 2650;

% Rotational speed of drill pipe
rpm = 150; % [50; 100; 150; 200];

% Cross fit, PAC4
k_Cr = 0.04119;
n_Cr =	0.4662;
ny_0 =	0.2058; 
ny_inf = 1.13e-08; 

% Rheological time scales
lambda_Maxwell_1 = 5; % [1; 10];
lambda_Maxwell_2 = 1.82e-3;
lambda_3ITT_1 = 15; % [8; 15; 20; 28]
lambda_3ITT_2 = 600;

% Colors
% http://se.mathworks.com/help/matlab/ref/colorspec.html
colorlist = {'b' ... % Blue [0 0 1]
    'g' ... % Green [0 1 1]
    'c' ... % Cyan [0 1 1]
    'm' ... % Magenta [1 0 1]
    'y'}; % Yellow [1 1 0]

% Particle diameter [mm] --> [m]
d_p = [0.01; 0.1; 1; 10] /1000;

% Settling distance [m]
H = max(d_h);

% Stokes settling velocity guess (m/s)
v_set = [4e-7; 4e-5; 4e-3; 7e-1];

%% Create figures

% Main Flow Scale - 3ITT
fig_MFS_3ITT  = figure('color','w');
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [8e-2 1e4],...
    'ylim', [1e-4 1e4],...
    'box','on');

xlabel('Re_G [-]');
ylabel('De_L_e_g_e_n_d [-]');


% Main Flow Scale - Maxwell
fig_MFS_Maxwell = figure('color','w');
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [8e-2 1e4],...
    'ylim', [1e-4 1e4],...
    'box','on');
    % ,...

xlabel('Re_G [-]');
ylabel('De_L_e_g_e_n_d [-]');

% Particle Scale - 3ITT
fig_PS_3ITT = figure('color','w');
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-8 1e2],...
    'ylim', [1e-4 1e4],...
    'box','on');
    % ,...

xlabel('Re_G [-]');
ylabel('De_L_e_g_e_n_d [-]');

% Particle Scale - Maxwell
fig_PS_Maxwell = figure('color','w');
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-8 1e2],...
    'ylim', [1e-4 1e4],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('Re_G [-]');
ylabel('De_L_e_g_e_n_d [-]');


%% Evaluate time scales & De

for i = 1:length(d_h)
    
    % Assign max. vol. flow rate for current d_h to vol. flow rate range
    VolFlowRate(4) = VolFlowRateMax(i);
    
    % Evaluate timescales & plot De vs. Re
    CompMainFlowScale(d_o(i), d_i(i), d_h(i), L, VolFlowRate, rpm);
end

CompParticleScale( d_p, H, v_set );

%% Format figures

% Main Flow Scale - 3ITT
figure(fig_MFS_3ITT);
legend('Mean flight time',...
    'Large Eddy Turnover time',...
    'Rotation',...
    'Shear Rate',...
    'Kolmogorov',...
    'location','northwest');

    % Plot De_crit
    hline(1e1,'k-');
    hline(1e-1,'k-');

    % h=fill([1e-2 1e-2 1e4 1e4], [0.1 10 10 0.1] ,'r');
    % set(h,'facealpha',.5);
    % set(h,'EdgeColor','none');
    

% Main Flow Scale - Maxwell
figure(fig_MFS_Maxwell);
legend('Mean flight time',...
    'Large Eddy Turnover time',...
    'Rotation',...
    'Shear Rate',...
    'Kolmogorov',...
    'location','northwest');

    % Plot De_crit
    hline(1e1,'k-');
    hline(1e-1,'k-');

    % h=fill([1e-2 1e-2 1e4 1e4], [0.1 10 10 0.1] ,'r');
    % set(h,'facealpha',.5);
    % set(h,'EdgeColor','none');

    
% Particle Scale - 3ITT
figure(fig_PS_3ITT);
legend('Settling time',...
    'Large Eddy Turnover time',...
    'Shear Rate',...
    'Stokes relaxation time',...
    'location','northwest');
    

    % Plot De_crit
    hline(1e1,'k-');
    hline(1e-1,'k-');
    
%     h=fill([1e-4 1e-4 1e4 1e4], [0.1 10 10 0.1] ,'r');
%     set(h,'facealpha',.5);
%     set(h,'EdgeColor','none');


% Particle Scale - Maxwell
figure(fig_PS_Maxwell);
legend('Settling time',...
    'Large Eddy Turnover time',...
    'Shear Rate',...
    'Stokes relaxation time',...
    'location','northwest');

    % Plot De_crit
    hline(1e1,'k-');
    hline(1e-1,'k-');
    
%     h=fill([1e-4 1e-4 1e4 1e4], [0.1 10 10 0.1] ,'r');
%     set(h,'facealpha',.5);
%     set(h,'EdgeColor','none');