%% Clean up
clear all;
close all;
clc;

col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 


addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';

%% Import SP data
% 
% [~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'SP_300_600','A3:F1724');
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
fig = figure('color','w','Units','centimeters','Position',[1 1 12 8.45]); % DIN A5 

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

% yyaxis left;
hold on;

% Set y-axis limits
ylim_min = 2e-2;
ylim_max = 3e-1;

% Format y-axis
ax1=gca;
set(gca,...
    'YScale','log',...
    'xlim', [0 1400],...
    'XTick', [0 200 400 600 800 1000 1200],...
    'ylim', [ylim_min ylim_max],...
    'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1],...
    'box','on');

xlabel('$t$ [s]','FontSize',10,'Interpreter','latex');
ylabel('$\eta$ [Pa.s]','FontSize',10,'Interpreter','latex');
grid('on');

plot([300 300],ylim,'k');
plot([900 900],ylim,'k');

% Expand Axes to Fill Figure (--> Minimum white space)
outerpos1 = ax1.OuterPosition; % [left bottom width height]


ti1 = ax1.TightInset; % [left bottom right top]


factor_hor = 0.01;
factor_ver = 0.01;

left = outerpos1(1) + ti1(1) + factor_hor;
bottom = outerpos1(2) + ti1(2) + factor_ver;

ax1_width = outerpos1(3) - ti1(1) - 2*factor_hor;
ax1_height = outerpos1(4) - ti1(2) - 2*factor_ver;

ax1.Position = [left bottom ax1_width ax1_height];



%% 300-600, PAC4

text(ax1,600,0.2,'PAC4','HorizontalAlignment','center','VerticalAlignment','Bottom','FontSize',14,'FontWeight','bold','Color',col_PAC4);

% Import UiS, PAC4 data

[~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'UiS_300_600','A3:F1256');

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

% Plot UiS data
% yyaxis left;
plot(ax1,Time, 0.1935*ones(length(Time)),'--','LineWidth',2,'Color',col_PAC4);
g(2) = plot(ax1,Time,Viscosity,'o','Color',col_PAC4,'markersize',MS,'MarkerFaceColor',MFC); % Viscosity
% yyaxis right;
% plot(Time,ShearStress,'o','markersize',MS,'MarkerFaceColor',MFC); % Shear stress


% Fit of growth function
Time_extract = Time(661:823)-Time(661);
Viscosity_extract = Viscosity(661:823);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(661);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_300_600_PAC4 = coefficients(3); 
TimeInterval = 0:0.5:150;
plot (ax1,t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_300_600_PAC4.*TimeInterval)),'-','Color','k','LineWidth',2);
txt=cat(2,'\lambda = ',num2str(round(lambda_3ITT_300_600_PAC4,1)),' s');
text(ax1,930,0.15,txt,'HorizontalAlignment','left','FontSize',10,'Color','k');


% Legend

% legend( g(1:2),...
%         'Lab1',...
%         'Lab2',...
%         'Location','southeast');
    
% Save timescales



%% 300-600, PAC2
text(ax1,600,0.051,'PAC2','HorizontalAlignment','center','VerticalAlignment','Bottom','FontSize',14,'FontWeight','bold','Color',col_PAC2);

[~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'PAC2_SH_300_600','B32:D1411');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Time = data(:,1);
ShearRate = data(:,2);
Viscosity = data(:,3);

% Clear temporary variables
clearvars data raw;

% Plot UiS data
% yyaxis left;
plot(ax1,Time,0.05*ones(length(Time)),'--','LineWidth',2,'Color',col_PAC2);
% g(3) = plot(ax1,Time,Viscosity,'o','Color',col_PAC2,'markersize',MS,'MarkerFaceColor','w'); % Viscosity
g(3) = plot(ax1,Time,[Viscosity(Time<=900); smooth(Viscosity(Time>900),9)],'o','Color',col_PAC2,'markersize',MS,'MarkerFaceColor','w'); % Viscosity
% yyaxis right;
% plot(Time,ShearStress,'o','markersize',MS,'MarkerFaceColor',MFC); % Shear stress


% % Fit of growth function
Time_extract = Time(890:1300)-Time(890);
Viscosity_extract = Viscosity(890:1300);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(890);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_300_600_PAC2 = coefficients(3); 
TimeInterval = 0:0.5:250;
plot (ax1,t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_300_600_PAC2.*TimeInterval)),'-','Color','k','LineWidth',2);
txt=cat(2,'\lambda = ',num2str(round(lambda_3ITT_300_600_PAC2,1)),' s');
text(ax1,1120,0.039,txt,'HorizontalAlignment','left','FontSize',10,'Color','k');



% Time(700)
% figure
% hold on;
% plot(Time(660:680),Viscosity(660:680))
% set(gca,'YScale','log');

Time_extract = Time(667:690)-Time(667);
Viscosity_extract = Viscosity(667:690);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'increasing');
coefficients = coeffvalues(fitresult);

t0      = Time(667);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_300_600_PAC2 = coefficients(3); 
TimeInterval = 0:0.5:13;

% Create smaller axes in bottom right and plot on it
ax3 = axes('Position',[.325 .52 .3 .3]); hold on; grid on;
plot(ax3,t0+Time_extract,smooth(Viscosity_extract,9),'o','Color',col_PAC2,'markersize',MS,'MarkerFaceColor','w');
plot (ax3,t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_300_600_PAC2.*TimeInterval)),'-','Color','k','LineWidth',2);
set(ax3,...
    'box','on',...
    'xlim',[t0 t0+TimeInterval(end)],...
    'ylim',[0.01 0.05],...
    'YScale','log',...
    'xtick',[],...
    'xticklabel',[],...
    'ytick',[],...
    'yticklabel',[],...
    'XColor','r',...
    'YColor','r');
% Add text with timescale
txt=cat(2,'\lambda = ',num2str(round(lambda_3ITT_300_600_PAC2,1)),' s');
text(ax3,908.5,0.02,txt,'HorizontalAlignment','left','FontSize',10,'Color','k');
rectangle(ax1,'Position',[900 0.0206 30 0.025],'EdgeColor','r','LineWidth',1);
xpos = [0.68 0.625];
ypos = [0.4 0.52];
annotation('line',xpos,ypos,'Color','r','LineWidth',1);
set(ax3,'linewidth',1)




%% Print as image

% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing


% Save Figure to File Format
path1 = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = '3ITT-300-600';
print(fig,[path1 fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[path1 fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 num2str(6)],'-dpng','-r600')
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file














%% Create figure
fig2 = figure('color','w','Units','centimeters','Position',[1 1 12 8.45]); % DIN A5 

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

% yyaxis left;
hold on;

% Set y-axis limits
ylim_min = 2e-2;
ylim_max = 3e-1;

% Format y-axis
ax1=gca;
set(gca,...
    'YScale','log',...
    'xlim', [0 650],...
    'XTick', [0 100 200 300 400 500 600],...
    'ylim', [ylim_min ylim_max],...
    'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1],...
    'box','on');

xlabel('$t$ [s]','FontSize',10,'Interpreter','latex');
ylabel('$\eta$ [Pa.s]','FontSize',10,'Interpreter','latex');
grid('on');


plot([100 100],[0.02 ylim_max],'k');
plot([150 150],[0.02 ylim_max],'k');

% Expand Axes to Fill Figure (--> Minimum white space)
outerpos1 = ax1.OuterPosition; % [left bottom width height]


ti1 = ax1.TightInset; % [left bottom right top]


factor_hor = 0.01;
factor_ver = 0.01;

left = outerpos1(1) + ti1(1) + factor_hor;
bottom = outerpos1(2) + ti1(2) + factor_ver;

ax1_width = outerpos1(3) - ti1(1) - 2*factor_hor;
ax1_height = outerpos1(4) - ti1(2) - 2*factor_ver;

ax1.Position = [left bottom ax1_width ax1_height];






%% 100-150, PAC4
text(ax1,600,0.21,'PAC4','HorizontalAlignment','center','VerticalAlignment','Bottom','FontSize',14,'FontWeight','bold','Color',col_PAC4);

% Import UiS, PAC4 data
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

% Plot UiS data
% yyaxis left;
plot(ax1,Time, 0.1935*ones(length(Time)),'--','LineWidth',2,'Color',col_PAC4); % Reference line
g(2) = plot(Time,Viscosity,'o','Color',col_PAC4,'markersize',MS,'MarkerFaceColor',MFC); % Viscosity
% yyaxis right;
% plot(Time,0.01981*ones(length(Time)),'k--','LineWidth',2); % Reference line
% plot(Time,ShearStress,'o','markersize',MS,'MarkerFaceColor',MFC); % Shear stress


% Fit of growth function
Time_extract = Time(662:785)-Time(662);
Viscosity_extract = Viscosity(662:785);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract,'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(662);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_100_50_PAC4     = coefficients(3); 
TimeInterval = 0:0.5:150;
plot (t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_100_50_PAC4.*TimeInterval)),'-','Color','k','LineWidth',2);
txt=cat(2,'\lambda = ',num2str(round(lambda_3ITT_100_50_PAC4,1)),' s');
text(180,0.26,txt,'HorizontalAlignment','left','FontSize',10,'Color','k');

%% 50-100, PAC2
text(ax1,600,0.052,'PAC2','HorizontalAlignment','center','VerticalAlignment','Bottom','FontSize',14,'FontWeight','bold','Color',col_PAC2);
[~, ~, raw] = xlsread([pwd '\3ITT_RRR.xlsx'],'PAC2_SH_100_50','B3:D822');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Time = data(:,1);
ShearRate = data(:,2);
Viscosity = data(:,3);

% Clear temporary variables
clearvars data raw;

% Plot UiS data
% yyaxis left;
plot(ax1,Time, 0.05*ones(length(Time)),'--','LineWidth',2,'Color',col_PAC2); % Reference 
g(4) = plot(ax1,Time,[smooth(Viscosity(Time<=100),9); Viscosity((Time>100)&(Time<=150)); smooth(Viscosity(Time>150),9)],'o','Color',col_PAC2,'markersize',MS,'MarkerFaceColor','w'); % Viscosity
% yyaxis right;
% plot(Time,ShearStress,'o','markersize',MS,'MarkerFaceColor',MFC); % Shear stress


% Fit of growth function
Time_extract = Time(111:250)-Time(111);
Viscosity_extract = Viscosity(111:250);
[fitresult, gof] = createFit(Time_extract, Viscosity_extract, 'increasing');
coefficients = coeffvalues(fitresult);

% e-function and time-constant
t0      = Time(111);
ny_0  = coefficients(1);
ny_inf  = coefficients(2);
lambda_3ITT_100_50_PAC2 = coefficients(3); 
TimeInterval = 0:0.5:150;
plot (ax1,t0+TimeInterval',ny_0+(ny_inf-ny_0).*(1-exp(-1./lambda_3ITT_100_50_PAC2.*TimeInterval)),'-','Color','k','LineWidth',2);
txt=cat(2,'\lambda = ',num2str(round(lambda_3ITT_100_50_PAC2,1)),' s');
text(ax1,180,0.071,txt,'HorizontalAlignment','left','FontSize',10,'Color','k');


% Save timescales

parentpath = cd(cd('..'));
path1 = [parentpath '\Timescales'];
save([path1 '\Timescales_3ITT.mat'],...
    'lambda_3ITT_300_600_PAC2',...
    'lambda_3ITT_300_600_PAC4',...
    'lambda_3ITT_100_50_PAC2',...
    'lambda_3ITT_100_50_PAC4');


% figure
% plot(Time(100:111),Viscosity(100:111))
% set(gca,'YScale','log');

    
%% Print as image

% Specify Figure Size and Page Size
set(fig2,'PaperPositionMode','auto') %set paper pos for printing


% Save Figure to File Format
path1 = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = '3ITT-100-50';
print(fig2,[path1 fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig2,[path1 fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig2,[path1 fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig2,[path1 num2str(7)],'-dpng','-r600')
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file