clear all;
close all;
clc;

addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';

col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 


%% Create figure

close all;
fig = figure('Units','centimeters','Position',[1 1 12 8.45]); % DIN A5 
hold on;

set(gcf,...
    'color','w');


%% PAC2

% Import SP data
[~, ~, raw] = xlsread( [pwd '\AS.xlsx'],'PAC2','A3:L28');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
MeasPts = data(:,1);
Time = data(:,2);
MeasPts1 = data(:,3);
Time1 = data(:,4);
AngularFrequency = data(:,5);
Strain = data(:,6);
ShearRate = data(:,7);
ShearStress = data(:,8);
StorageModulus = data(:,9);
LossModulus = data(:,10);
ComplexModulus = data(:,11);
ComplexViscosity = data(:,12);

% Clear temporary variables
clearvars data raw;

% Plot

MS = 6; % MarkerSize
MFC = 'k'; % MarkerFaceColor

plot(Strain,LossModulus,...
    '-^','markersize',MS,'color',col_PAC2,'MarkerFaceColor',col_PAC2);
plot(Strain,StorageModulus,...
    '-s','markersize',MS,'color',col_PAC2,'MarkerFaceColor',col_PAC2);


% Phase shift angle on second y-axis
yyaxis right;
plot(Strain,atan(LossModulus./StorageModulus).*180./pi,...
    '-','Color',col_PAC2,'linewidth',2);
set(gca,'Ycolor','k'); 
ylabel('$\phi$ [^o]','FontSize',10,'Interpreter','latex');
yyaxis left;



% % Sigve 2018 NTNU Petroleum lab
% [~, ~, raw] = xlsread( [pwd '\180205_PAC2_AS_Sigve.xlsx'],'1rad','AK36:AU55');
% 
% % Create output variable
% data = reshape([raw{:}],size(raw));
% 
% % Allocate imported array to column variable names
% MeasPts = data(:,1);
% Time = data(:,2);
% AngularFrequency = data(:,3);
% Strain = data(:,4);
% ShearRate = data(:,5);
% ShearStress = data(:,6);
% StorageModulus = data(:,7);
% LossModulus = data(:,8);
% ComplexModulus = data(:,9);
% ComplexViscosity = data(:,10);
% DeflAngle = data(:,11);
% 
% % Clear temporary variables
% clearvars data raw;
% 
% plot(Strain,LossModulus,...
%     '--r^','markersize',MS);
% plot(Strain,StorageModulus,...
%     '--rs','markersize',MS);


%% PAC4

% Import SP data
[~, ~, raw] = xlsread( [pwd '\AS.xlsx'],'PAC4_MK','A3:L28');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
MeasPts = data(:,1);
Time = data(:,2);
MeasPts1 = data(:,3);
Time1 = data(:,4);
AngularFrequency = data(:,5);
Strain = data(:,6);
ShearRate = data(:,7);
ShearStress = data(:,8);
StorageModulus = data(:,9);
LossModulus = data(:,10);
ComplexModulus = data(:,11);
ComplexViscosity = data(:,12);

% Clear temporary variables
clearvars data raw;



% Plot

MS = 6; % MarkerSize


plot(Strain,LossModulus,...
    '-^','markersize',MS,'MarkerFaceColor',col_PAC4,'color',col_PAC4);
plot(Strain,StorageModulus,...
    '-s','markersize',MS,'MarkerFaceColor',col_PAC4,'color',col_PAC4);


% Phase shift angle on second y-axis
yyaxis right;
plot(Strain,atan(LossModulus./StorageModulus).*180./pi,...
    '-','Color',col_PAC4,'linewidth',2);
set(gca,'Ycolor','k'); 
ylabel('$\phi$ [$^o$]','FontSize',10,'Interpreter','latex');
yyaxis left;



%% Formating

xlabel('$\gamma$ [\%]','FontSize',10,'Interpreter','latex');
ylabel('$G''$, $G''''$ [Pa]','FontSize',10,'Interpreter','latex');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-1 1e3],...
    'ylim', [1e-3 2e-1],...
    'XTick', [1e-2  1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [1e-3 1e-2 1e-1 1e-0],...
    'box','on');


%% Legend
str={'G'''', PAC2', ...    
    'G'', PAC2',...    
    'G'''', PAC4', ...    
    'G'', PAC4',...
    '\phi, PAC2',...
    '\phi, PAC4'}';
    
% h_legend = legend(...
%     'G'''', PAC2', ...    
%     'G'', PAC2',...    
%     'G'''', PAC4', ...    
%     'G'', PAC4',...
%     '\phi, PAC2',...
%     '\phi, PAC4');
% set(h_legend,'Location','southeast','FontSize',10);
% legend boxoff;
columnlegend(3, str,'location','southeast');
pos = get(legend, 'position');
set(legend, 'position', [1.05*pos(1) -0.9*pos(2) pos(3) 1.4*pos(1)]);


%% Additional information

% linecolor = [161 149 137]; %SINTEF grey
% addpath('M:\data\SW working folders\MATLAB\m-files');
% vline([100 1000],{'k--' 'k--'},{'$\phi = 73°$' '$\phi = 86°$'});
% hline([LossModulus(21) LossModulus(26) StorageModulus(21) StorageModulus(26)],{'k--' 'k--' 'k--' 'k--'},{'' '' '' ''});


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

%fig.PaperPosition = [1 1 12 8.45]

%fig_pos = fig.PaperPosition;
%fig.PaperUnits = 'centimeters';
%fig.PaperType = '<custom>';
%fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.PaperOrientation = 'landscape';
fig.Position
fig.PaperPosition
fig.OuterPosition
fig.PaperSize

% Save Figure to File Format
fig_path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = 'AS';
print(fig,[fig_path fig_name],'-dpdf');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig,[fig_path num2str(3)],'-dpng','-r600') 
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file