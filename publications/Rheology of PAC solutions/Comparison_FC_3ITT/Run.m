%% Clean up
clear all;
close all;
clc;



col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 



%% Create figure
fig = figure('color','w','Units','centimeters','Position',[1 1 12 8.45]); % DIN A5 

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

% Format of figure

% yyaxis left;
hold on;

% Set y-axis limits
ylim_min = 0.8;
ylim_max = 1.2;

% Format y-axis
set(gca,...
    'XScale','lin',...
    'YScale','lin',...
    'xlim', [0 1400],...
    'XTick', [0 200 400 600 800 1000 1200],...
    'YTick', [0.8 0.9 1.0 1.1 1.2],...
    'ylim', [ylim_min ylim_max],...
    'box','on');

xlabel('$t$ (300-600-500 test) [s]','FontSize',10,'Interpreter','latex');
ylabel('$\eta_{3ITT}/\eta_{FC}$ [-]','FontSize',10,'Interpreter','latex');
grid('on');


plot([300 300],[ylim_min 1.1],'k');
plot([900 900],ylim,'k');

% Second axis
ax1 = gca; ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'XScale','lin',...
    'YScale','lin',...
    'xlim', [-750 650],...
    'ylim', ylim,...
    'FontSize',10,...
    'XTick',[0 150 250 450],...
    'YTick', []);

% ylabel('Relaxation time [s]');

% Expand Axes to Fill Figure (--> Minimum white space)
outerpos1 = ax1.OuterPosition; % [left bottom width height]
outerpos2 = ax2.OuterPosition; % [left bottom width height]

ti1 = ax1.TightInset; % [left bottom right top]
ti2 = ax2.TightInset; % [left bottom right top]

factor_hor = 0.01;
factor_ver = 0.01;

left = outerpos1(1) + ti1(1) + factor_hor;
bottom = outerpos1(2) + ti1(2) + factor_ver;

ax1_width = outerpos1(3) - ti1(1) - ti2(3) - 2*factor_hor;
ax1_height = outerpos1(4) - ti1(2) - ti2(4) - 2*factor_ver;

ax1.Position = [left bottom ax1_width ax1_height];
ax2.Position = [left bottom ax1_width ax1_height];

hold on;

% Origin for second axis
plot(ax2,[0 0],[0.95 max(ylim)],'k-','LineWidth',0.5);
plot(ax2,[100 100],[0.95 max(ylim)],'k-','LineWidth',0.5);
% second x-axis title
xlabel('$t$ (100-50-500 test) [s]','FontSize',10,'Interpreter','latex');
vec_pos = get(get(ax2, 'XLabel'), 'Position');
set(get(ax2, 'XLabel'), 'Position', vec_pos + [-220 -0.03 0]);

ax1.Layer = 'top';
ax2.Layer = 'top';

%% UIS 300-600 PAC4

parentpath = cd(cd('..'));
path1 = [parentpath '\3ITT-RRR\']; 
addpath(path1);

%Import UiS 3ITT data
[~, ~, raw] = xlsread([path1 '3ITT_RRR.xlsx'],'UiS_300_600','A3:F1111');

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

% Import UiS FC data
path2 = [parentpath '\FC\Data\Lab2_UiS\'];
[~, ~, raw] = xlsread([path2 '160629_PAC4_UiS.xlsx'],'A3:B28');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
UiS_FC_ShearRate = data(:,1);
UiS_FC_Viscosity = data(:,2);

% Clear temporary variables
clearvars data raw;

% Plot Ratio of 3ITT and FC viscositities using interpolation or raw data
plot(ax1,Time, Viscosity./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, ShearRate,'cubic','extrap'),'o','Color',col_PAC4);

% Fit of growth function - High shear rate interval
Time_extract = Time(66:660)-Time(66);
Viscosity_extract = Viscosity(66:660);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'decreasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(65);
ny_inf  = coefficients(1);
ny_0  = coefficients(2);
lambda_3ITT_300_600_PAC4_HighShearRate     = coefficients(3); 
TimeInterval = 0:0.5:600;
plot (ax1,t0+TimeInterval',(ny_inf+(ny_0-ny_inf).*exp(-1./lambda_3ITT_300_600_PAC4_HighShearRate.*TimeInterval))./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, mean(ShearRate(66:660))),'k-','LineWidth',2);
txt=cat(2,'\lambda_3_0_0_-_6_0_0_-_._._. = ',num2str(round(lambda_3ITT_300_600_PAC4_HighShearRate,0)),' s');
text(ax1,620,0.905,txt,'HorizontalAlignment','center','VerticalAlignment','Top','FontSize',10,'Color','black');

% Fit of growth function - Low shear rate interval
Time_extract = Time(720:1109)-Time(720);
Viscosity_extract = Viscosity(720:1109);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(719);
ny_inf  = coefficients(2);
ny_0  = coefficients(1);
lambda_3ITT_300_600_PAC4_LowShearRate     = coefficients(3); 
TimeInterval = 0:0.5:600;
plot (ax1,t0+TimeInterval',(ny_inf+(ny_0-ny_inf).*exp(-1./lambda_3ITT_300_600_PAC4_LowShearRate.*TimeInterval))./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, mean(ShearRate(720:1109))),'k-','LineWidth',2);
txt=cat(2,'\lambda_3_0_0_-_6_0_0_-_._._. = ',num2str(round(lambda_3ITT_300_600_PAC4_LowShearRate,0)),' s');
text(ax1,1390,0.88,txt,'HorizontalAlignment','right','FontSize',10,'Color','black');


%% UIS 100-50 PAC4

%Import UiS 3ITT data
[~, ~, raw] = xlsread([path1 '3ITT_RRR.xlsx'],'UiS_100_50','A3:F1111');

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

% % Import UiS FC data
% path = [parentpath '\FC\Data\Lab2_UiS\'];
% [~, ~, raw] = xlsread([path '160629_PAC4_UiS.xlsx'],'A3:B28');
% 
% % Create output variable
% data = reshape([raw{:}],size(raw));
% 
% % Allocate imported array to column variable names
% UiS_FC_ShearRate = data(:,1);
% UiS_FC_Viscosity = data(:,2);
% 
% % Clear temporary variables
% clearvars data raw;

% Plot Ratio of 3ITT and FC viscositities using interpolation or raw data
plot(ax2,Time, Viscosity./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, ShearRate,'cubic','extrap'),'o','Color',col_PAC4);

% Fit of growth function - Low shear rate interval
Time_extract = Time(829:1109)-Time(829);
Viscosity_extract = Viscosity(829:1109);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'decreasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(829);
ny_inf  = coefficients(1);
ny_0  = coefficients(2);
lambda_3ITT_100_50_PAC4_LowShearRate = coefficients(3); 
TimeInterval = 0:0.5:300;
plot (ax2,t0+TimeInterval',(ny_inf+(ny_0-ny_inf).*exp(-1./lambda_3ITT_100_50_PAC4_LowShearRate.*TimeInterval))./interp1(UiS_FC_ShearRate,UiS_FC_Viscosity, mean(ShearRate(829:1109))),'k-','LineWidth',2);
txt=cat(2,'\lambda_1_0_0_-_5_0_-_._._. = ',num2str(round(lambda_3ITT_100_50_PAC4_LowShearRate,0)),' s');
text(ax2,400,1.15,txt,'HorizontalAlignment','center','FontSize',10,'Color','black');



%% Plot references
Time = get(ax1,'xlim');
plot(ax1,Time,ones(length(Time)),'--','Color',col_PAC4,'LineWidth',2);


    
%% Save timescales
parentpath = cd(cd('..'));
path1 = [parentpath '\Timescales'];
save([path1 '\Timescales_3ITT_Large.mat'],... 
    'lambda_3ITT_300_600_PAC4_HighShearRate',...
    'lambda_3ITT_300_600_PAC4_LowShearRate',...
    'lambda_3ITT_100_50_PAC4_LowShearRate');
    

%% Print as image

% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperUnits = 'centimeters';
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path1 = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = '3ITT_vs_FC';
print(fig,[path1 fig_name],'-dpdf');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[path1 fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 num2str(8)],'-dpng','-r600');
% saveas(fig, [path fig_name],'meta') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file