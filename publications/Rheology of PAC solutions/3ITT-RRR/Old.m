%% Clean up
clear all;
close all;
clc;

%% Import SP data
% [~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'SP_100_50','A3:F1181');
% 
% % Create output variable
% data = reshape([raw{:}],size(raw));
% 
% % Allocate imported array to column variable names
% Time = data(:,1);
% ShearRate = data(:,2);
% ShearStress = data(:,3);
% Viscosity = data(:,4);
% Speed = data(:,5);
% Torque = data(:,6);
% 
% % Clear temporary variables
% clearvars data raw;

%% Create figure
fig = figure('color','w','Units','centimeters','Position',[1 1 21 14.8]); % DIN A5 

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

%% Format of figure

% yyaxis left;
hold on;

% Set y-axis limits
ylim_min = 2e-2;
ylim_max = 3e-1;

% Format y-axis
set(gca,...
    'YScale','log',...
    'xlim', [0 863],...
    'ylim', [ylim_min ylim_max],...
    'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1],...
    'box','on',...
    'XAxisLocation','top');

xlabel('Time [s]');
grid('on');

% Interval 2
HSR = area([100 150], [ylim_max ylim_max]);
HSR.FaceColor = [0.8 0.8 0.8];
HSR.EdgeColor = 'none';
HSR.FaceAlpha = 0.5;

% Interval shear rates labeling
y1 = 0.055;
dy1=0.03;
dy2=0.01;

x1 = 50;
txt1 = '($\dot{\gamma}~=~0.1~s^{-1}$)';
text(x1,y1,txt1,'HorizontalAlignment','center','FontSize',10,'interpreter','latex');
text(x1,y1+dy2,'Rotation','HorizontalAlignment','center','FontSize',10);

x1 = 125;
txt1 = '($\dot{\gamma}~=~1200~s^{-1}$)';
text(x1,y1+dy1,txt1,'HorizontalAlignment','center','FontSize',10,'interpreter','latex');
text(x1,y1+dy1+dy2,'Rotation','HorizontalAlignment','center','FontSize',10);

x1 = 350;
txt1 = '($\dot{\gamma}~=~0.1~s^{-1})$';
text(x1,y1,txt1,'HorizontalAlignment','center','FontSize',10,'interpreter','latex');
text(x1,y1+dy2,'Rotation','HorizontalAlignment','center','FontSize',10);


%% Viscosity as function of time
% g(1) = plot(Time,Viscosity,'kd','markersize',MS,'MarkerFaceColor',MFC); % Viscosity
% 
% % Fit of growth function
% Time_extract = Time(128:419)-Time(128);
% Viscosity_extract = Viscosity(128:419);
% [fitresult, gof] = createFit(Time_extract, Viscosity_extract,'increasing');
% coefficients = coeffvalues(fitresult);
% 
% % e-function and time-constant
% t0      = Time(128);
% ny_0  = coefficients(1);
% ny_inf  = coefficients(2);
% lambda_3ITT_100_50_Lab1     = coefficients(3); 
% TimeInterval = 0:0.5:150;
% plot (t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_100_50_Lab1.*TimeInterval)),'b-','LineWidth',2);
% txt=cat(2,'\lambda_L_a_b_1 = ',num2str(round(lambda_3ITT_100_50_Lab1,1)),' s');
% text(230,0.27,txt,'HorizontalAlignment','center','FontSize',10,'Color','blue');

% Formating
ylabel('Apparent Viscosity [Pa·s]');

%% Shear stress as function of time

% yyaxis right;
% plot(Time,ShearStress,'d','markersize',MS,'MarkerFaceColor',MFC); % Shear stress
% 
% % Formating
% ylabel('Stress [Pa]');
% 
% % Set y-axis limits
% ylim_min = 1e-2;
% ylim_max = 4e1;
% 
% % Format y-axis
% set(gca,...
%     'YScale','log',...
%     'ylim', [ylim_min ylim_max],...
%     'YTick', [1e-2 1e-1 1e-0 1e1 4e1],...
%     'FontSize',24);
% 


%% Import UiS data
[~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'UiS_100_50','A3:F1256');

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



%% Plot UiS data
% yyaxis left;
plot(Time,0.2*ones(length(Time)),'k--','LineWidth',2); % Reference line
g(2) = plot(750+Time,Viscosity,'go','markersize',MS,'MarkerFaceColor',MFC); % Viscosity
% yyaxis right;
% plot(Time,0.01981*ones(length(Time)),'k--','LineWidth',2); % Reference line
% plot(Time,ShearStress,'o','markersize',MS,'MarkerFaceColor',MFC); % Shear stress


% Fit of growth function
Time_extract = Time(662:785)-Time(662);
Viscosity_extract = Viscosity(662:785);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract,'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = 750+Time(662);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_100_50_Lab2     = coefficients(3); 
TimeInterval = 0:0.5:150;
plot (t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_100_50_Lab2.*TimeInterval)),'b-','LineWidth',2);
txt=cat(2,'\lambda_L_a_b_2 = ',num2str(round(lambda_3ITT_100_50_Lab2,1)),' s');
text(235,0.18,txt,'HorizontalAlignment','center','FontSize',10,'Color','blue');


%% Legend

% legend( g(1:2),...
%         'Lab1',...
%         'Lab2',...
%         'Location','southeast');
    
%% Save timescales
parentpath = cd(cd('..'));
path = [parentpath '\Timescales'];
save([path '\Timescales_3ITT_100_50.mat'],...
    'lambda_3ITT_100_50_Lab2');

    
    


%% Print as image

% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.01;
factor_ver = 0.01;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];




% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub1\Figures\';
fig_name = '3ITT-RRR_100_50';
print(fig,[path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file
