% Comparison of 3ITT-RRR & FC

clear all;
close all;
clc;

%% Create figure
fig = figure('color','w'); % ('units','normalized','outerposition',[0 0 1 1]);
hold on;

path = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\MATLAB\AdWell\Rheology\PAC\3ITT-RRR\';

%% Formating

xlabel('Time');
ylabel('Apparent viscosity ratio \eta_3_I_T_T / \eta_F_C');
grid('on');

ylim_min = 0.8;
ylim_max = 1.4;

set(gca,...
    'xlim', [0 863],...
    'ylim', [ylim_min ylim_max],...
    'box','on',...
    'FontSize',24);


%     'XScale','log',...
%     'YScale','log',...

%     'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
%     'YTick', [2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1],...

% xgrid = get(gca,'XGridHandle');  % or: hAxes.XGridHandle
% ygrid = get(gca,'YGridHandle');  % or: hAxes.XGridHandle
% set(xgrid,'Color','k', 'GridLineStyle','-', 'LineWidth',1, 'MajorMinor','majorandminor');
% set(ygrid,'Color','k', 'GridLineStyle','-', 'LineWidth',1, 'MajorMinor','majorandminor');

HSR = area([100 150], [ylim_max ylim_max]);
HSR.FaceColor = [0.8 0.8 0.8];
HSR.EdgeColor = 'none';
HSR.FaceAlpha = 0.5;


%% Import SP 3ITT 100_50 data
[~, ~, raw] = xlsread([path '3ITT_RRR.xlsx'],'SP_100_50','A3:F1181');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Time = data(:,1);
ShearRate = data(:,2);
ShearStress = data(:,3);
Viscosity = data(:,4);
Speed = data(:,5);
Torque = data(:,6);

% Clear temporary variables
clearvars data raw;


%% Plot Ratio of 3ITT and FC viscositities using interpolation or raw data

% Lab 2 FC data

load FC_Lab1;
% g(1) = plot(Time, Viscosity./interp1(Up_SR,Up_AV_M, ShearRate),'k--');
g(2) = plot(Time, Viscosity./interp1(Up_SR,mean([Up_AV_M flip(Down_AV_M)],2), ShearRate),'kd');
% g(3) = plot(Time, Viscosity./interp1(Down_SR,Down_AV_M, ShearRate),'k:');

% addpath 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\MATLAB\AdWell\Rheology\PAC\FC\Data\Lab2_UiS'
% load Lab2_Scatter2016;
% g(4) = plot(Time, Viscosity./interp1(Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2), ShearRate,'cubic','extrap'),'r--');

% help = Viscosity./interp1(Up_SR,mean([Up_AV_M flip(Down_AV_M)],2), ShearRate);
% 
% 
% % Fit of growth function - High shear rate interval
% Time_extract = Time(31:119)-Time(31);
% Viscosity_extract = help(31:119);
% [fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'decreasing');
% coefficients = coeffvalues(fitresult);
% 
% % e-function and time-constant
% t0      = Time(30);
% ny_inf  = coefficients(1);
% ny_0  = coefficients(2);
% tau     = coefficients(3); 
% TimeInterval = 0:0.5:50;
% % 
% % figure
% % hold on
% % scatter (t0+Time_extract, Viscosity_extract)
% 
% plot (t0+TimeInterval',(ny_inf+(ny_0-ny_inf).*exp(-1./tau.*TimeInterval))./interp1(Up_SR,Up_AV_M, mean(ShearRate(66:660))),'r-','LineWidth',2);
% txt=cat(2,'\tau_L_a_b_1 = ',num2str(tau),' s');
% text(150,1.02,txt,'HorizontalAlignment','center','FontSize',18,'Color','red');


%% Import UiS 3ITT data
[~, ~, raw] = xlsread([path '3ITT_RRR.xlsx'],'UiS_100_50','A3:F1256');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Time = data(:,1);
ShearRate = data(:,2);
ShearStress = data(:,3);
Viscosity = data(:,4);
Speed = data(:,5);
Torque = data(:,6);

% Clear temporary variables
clearvars data raw;

%% Import UiS FC data
path = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\MATLAB\AdWell\Rheology\PAC\FC\Data\Lab2_UiS\';
[~, ~, raw] = xlsread([path '160629_PAC4_UiS.xlsx'],'A3:B28');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
UiS_FC_ShearRate = data(:,1);
UiS_FC_Viscosity = data(:,2);

% Clear temporary variables
clearvars data raw;

% help = Viscosity./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, ShearRate,'cubic','extrap');
% 
% 
% % Fit of growth function - High shear rate interval
% Time_extract = Time(100:660)-Time(68);
% Viscosity_extract = help(100:660);
% [fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'decreasing');
% coefficients = coeffvalues(fitresult);
% 
% 
% % e-function and time-constant
% t0      = Time(67);
% ny_inf  = coefficients(1);
% ny_0  = coefficients(2);
% tau     = coefficients(3); 
% TimeInterval = 0:0.5:50;
% 
% % figure
% % hold on
% % scatter (t0+Time_extract, Viscosity_extract)
% 
% 
% plot (t0+TimeInterval',(ny_inf+(ny_0-ny_inf).*exp(-1./tau.*TimeInterval))./interp1(Up_SR,Up_AV_M, mean(ShearRate(66:660))),'r-','LineWidth',2);
% txt=cat(2,'\tau_L_a_b_1 = ',num2str(tau),' s');
% text(150,0.92,txt,'HorizontalAlignment','center','FontSize',18,'Color','red');





%% Plot Ratio of 3ITT and FC viscositities using interpolation or raw data

g(4) = plot(Time, Viscosity./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, ShearRate,'cubic','extrap'),'ko');

%% Plot references & intervalls

plot(Time,ones(length(Time)),'k--','LineWidth',2);

x1 = 50;
y1 = 0.85; 
txt1 = '$\dot{\gamma}~=~0.1~s^{-1}$';
text(x1,y1,txt1,'HorizontalAlignment','center','FontSize',18,'interpreter','latex');

x1 = 125;
y1 = 0.85;
txt1 = '$\dot{\gamma}~=~1200~s^{-1}$';
text(x1,y1,txt1,'HorizontalAlignment','center','FontSize',18,'interpreter','latex');

x1 = 275;
y1 = 0.85;
txt1 = '$\dot{\gamma}~=~0.1~s^{-1}$';
text(x1,y1,txt1,'HorizontalAlignment','center','FontSize',18,'interpreter','latex');

legend( g(1:4),...
        '\eta_F_C = Mean of all upward sweeps (Lab1)',...
        '\eta_F_C = Total mean of all upward and downward sweeps (Lab1)',...
        '\eta_F_C = Mean of all downward sweeps (Lab1)',...
        '\eta_F_C = Total mean of all upward and downward sweeps (Lab2)',...
        'Location','southeast');
    
    
%% Print as image

% Adjust size to screen size
screen_size = get(0, 'ScreenSize');
origSize = get(fig, 'Position'); % grab original on screen size
offset_hor = 300;
offset_ver = 80;
set(fig, 'Position', [offset_hor 0 screen_size(3)-2*offset_hor screen_size(4)-offset_ver ] ); %set to scren size


% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.003;
factor_ver = 0.01;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];

% Re-Adjust size
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
% set(fig,'Units','centimeters','Position',[0 0 32 18])
% set(fig,'Position', origSize) %set back to original dimensions

% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path = 'M:\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper\Figures\';
fig_name = 'Figure7a_3ITT_vs_Flowcurves';
print(fig,[path fig_name],'-dpdf');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig,[path fig_nam],'-dpng');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file