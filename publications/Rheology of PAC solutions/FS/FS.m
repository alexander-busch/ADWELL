clear all;
close all;
clc;

addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';
path1 = 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\PAC\Resub1\FS\';

col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 


%% Create figure

close all;
fig = figure('color','w','Units','centimeters','Position',[1 1 12 8.45]);

hold on;

% Formating

xlabel('$\omega$ [rad/s]','Interpreter','latex');
ylabel('$G''$, $G''''$ [Pa], $\lambda_{Max}$ [s]','Interpreter','latex');
grid('on');
box on;

ax1 = gca;  

xmin = 1e-0;
xmax = 1e+2;
ymin = 6e-3;
ymax = 2e+1;

set(ax1,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [xmin xmax],...
    'ylim', [ymin ymax]);

% Second axis
% ax2 = axes('Position',get(ax1,'Position'),...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'XColor','k','YColor','k',...
%     'XScale','log',...
%     'YScale','log',...
%     'xlim', [xmin xmax],...
%     'ylim', [ymin ymax],...
%     'FontSize',10);
% xlabel('Shear rate [1/s]');
% ylabel('Relaxation time \lambda_M_a_x [s]');

hold on;

%% PAC2 

% Import data
% [~, ~, raw] = xlsread([path 'FS.xlsx'],'PAC2','A4:I24'); % Sigve
[~, ~, raw] = xlsread([path1 'FS.xlsx'],'PAC2','A35:I45'); % Milad


% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
MeasPts = data(:,1);
AngularFrequency = data(:,2);
StorageModulus = data(:,3);
LossModulus = data(:,4);
ComplexModulus = data(:,5);
DynamicViscosity = data(:,6);
ComplexViscosity = data(:,7);
MaxwellRelaxationTime_PAC2 = data(:,8);
ShearRate = data(:,9);

% Clear temporary variables
clearvars data raw;

% Plot

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

g(1) = plot(ax1,AngularFrequency,LossModulus,...
    '-^','Color',col_PAC2,'markersize',MS,'MarkerFaceColor',col_PAC2);
g(2) = plot(ax1,AngularFrequency,StorageModulus,...
    '-s','Color',col_PAC2,'markersize',MS,'MarkerFaceColor',col_PAC2);
g(3) = plot(ax1,AngularFrequency,MaxwellRelaxationTime_PAC2,...
    '-','Color',col_PAC2,'markersize',MS,'MarkerFaceColor',col_PAC2,'LineWidth',2);




AngularFrequency_PAC2=AngularFrequency;
LossModulus_PAC2=LossModulus;
StorageModulus_PAC2=StorageModulus;
ComplexViscosity_PAC2=ComplexViscosity;

% Export Angular Frequency & Complex Viscosity for Cox-Merz 
parentpath = cd(cd('..'));
path2 = [parentpath '\NormalStressDifferences'];
save([path2 '\FS_PAC2.mat'],...
    'AngularFrequency_PAC2',...
	'LossModulus_PAC2',...
    'StorageModulus_PAC2',...
    'ComplexViscosity_PAC2');


% % Sigve 2018 NTNU Petroleum lab
% [~, ~, raw] = xlsread( [pwd '\180205_PAC2_FS_Sigve.xlsx'],'Ark1','I30:O40');
% % [~, ~, raw] = xlsread( [pwd '\180205_PAC2_FS_Sigve.xlsx'],'Ark1','A30:G40');
% 
% % Create output variable
% data = reshape([raw{:}],size(raw));
% 
% % Allocate imported array to column variable names
% MeasPts = data(:,1);
% AngularFrequency = data(:,2);
% StorageModulus = data(:,3);
% LossModulus = data(:,4);
% ComplexModulus = data(:,5);
% ComplexViscosity = data(:,6);
% Strain = data(:,7);
% 
% % Clear temporary variables
% clearvars data raw;
% 
% plot(AngularFrequency,LossModulus,...
%     '--b^','markersize',MS);
% plot(AngularFrequency,StorageModulus,...
%     '--bs','markersize',MS);

%% PAC 4 (SIntef Petroleum)

% % Import SP data
% [~, ~, raw] = xlsread([path 'FS.xlsx'],'SP','A3:I21');
% 
% % Create output variable
% data = reshape([raw{:}],size(raw));
% 
% % Allocate imported array to column variable names
% MeasPts = data(:,1);
% AngularFrequency = data(:,2);
% StorageModulus = data(:,3);
% LossModulus = data(:,4);
% ComplexModulus = data(:,5);
% DynamicViscosity = data(:,6);
% ComplexViscosity = data(:,7);
% MaxwellRelaxationTime = data(:,8);
% ShearRate = data(:,9);
% 
% % Clear temporary variables
% clearvars data raw;
% 
% % Plot 
% 
% MS = 6; % MarkerSize
% MFC = 'k'; % MarkerFaceColor
% 
% g(1) = plot(AngularFrequency,LossModulus,...
%     ':g^','markersize',MS,'MarkerFaceColor',MFC);
% g(2) = plot(AngularFrequency,StorageModulus,...
%     ':gs','markersize',MS,'MarkerFaceColor',MFC);
% g(3) = plot(AngularFrequency,MaxwellRelaxationTime,...
%     '-','Color','blue','markersize',MS,'MarkerFaceColor',MFC);

%% PAC 4 (UiS)

% Import UiS data
% [~, ~, raw] = xlsread([path 'FS.xlsx'],'UiS','A3:I23'); % Strain 0.2 %
% [~, ~, raw] = xlsread([path 'FS.xlsx'],'PAC4','A3:I13'); % Strain 1%
[~, ~, raw] = xlsread([path1 'FS.xlsx'],'PAC4','A18:I28'); % Strain 0.2 %


% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
MeasPts = data(:,1);
AngularFrequency = data(:,2);
StorageModulus = data(:,3);
LossModulus = data(:,4);
ComplexModulus = data(:,5);
DynamicViscosity = data(:,6);
ComplexViscosity = data(:,7);
MaxwellRelaxationTime_PAC4 = data(:,8);
ShearRate = data(:,9);

% Clear temporary variables
clearvars data raw;


% Plot

MS = 6; % MarkerSize

g(4) = plot(ax1,AngularFrequency,LossModulus,...
    '-^','Color',col_PAC4,'markersize',MS,'MarkerFaceColor',col_PAC4);
g(5) = plot(ax1,AngularFrequency,StorageModulus,...
    '-s','Color',col_PAC4,'markersize',MS,'MarkerFaceColor',col_PAC4);
g(6) = plot(ax1,AngularFrequency,MaxwellRelaxationTime_PAC4,...
    '-','Color',col_PAC4,'markersize',MS,'MarkerFaceColor',col_PAC4,'LineWidth',2);


AngularFrequency_PAC4=AngularFrequency;
LossModulus_PAC4=LossModulus;
StorageModulus_PAC4=StorageModulus;
ComplexViscosity_PAC4=ComplexViscosity;

% Export Angular Frequency & Complex Viscosity for Cox-Merz 
parentpath = cd(cd('..'));
path2 = [parentpath '\NormalStressDifferences'];
save([path2 '\FS_PAC4.mat'],...
    'AngularFrequency_PAC4',...
	'LossModulus_PAC4',...
    'StorageModulus_PAC4',...
    'ComplexViscosity_PAC4');




%% Legend
% h_legend = legend( g(1:6),...
%     'G'''', PAC2', ...    
%     'G'', PAC2',...    
%     '\lambda_M_a_x, PAC2',...   
%     'G'''', PAC2', ...    
%     'G'', PAC2',...    
%     '\lambda_M_a_x, PAC2');
% set(h_legend,'Location','northwest','FontSize',10);


%leg_pos=get(h_legend,'Position');
%set(h_legend,'Position',[0.43 0.11 leg_pos(3) leg_pos(4)],'FontSize',10);
%get(h_legend,'fontname');
str =  {'G'''', PAC2', ...    
    'G'', PAC2',...    
    '\lambda_M_a_x, PAC2',...   
    'G'''', PAC2', ...    
    'G'', PAC2',...    
    '\lambda_M_a_x, PAC2'}';
%columnlegend(3, str,'location','northwest');


%% Time scales

addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\viscosity-based'

% Power law fit, does not really fit...
MRT = createFit_PowerLaw(AngularFrequency,MaxwellRelaxationTime_PAC4);
close gcf;

parentpath = cd(cd('..'));
path1 = [parentpath '\Timescales'];
save([path1 '\Timescales_FS.mat'],...
    'AngularFrequency_PAC2',...
    'MaxwellRelaxationTime_PAC2',...
    'AngularFrequency_PAC4',...
    'MaxwellRelaxationTime_PAC4');


%% Print as image

        
% Expand Axes to Fill Figure (--> Minimum white space)
outerpos1 = ax1.OuterPosition; % [left bottom width height]
%outerpos2 = ax2.OuterPosition; % [left bottom width height]

ti1 = ax1.TightInset; % [left bottom right top]
%ti2 = ax2.TightInset; % [left bottom right top]

factor_hor = 0.01;
factor_ver = 0.01;

left = outerpos1(1) + ti1(1) + factor_hor;
bottom = outerpos1(2) + ti1(2) + factor_ver;

ax1_width = outerpos1(3) - ti1(1) - 2*factor_hor; %- ti2(3)
ax1_height = outerpos1(4) - ti1(2) - 2*factor_ver; %- ti2(4) 

ax1.Position = [left bottom ax1_width ax1_height];
% ax2.Position = [left bottom ax1_width ax1_height];

% Re-Adjust size
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
% set(fig,'Units','centimeters','Position',[0 0 32 18])
% set(fig,'Position', origSize) %set back to original dimensions
columnlegend(3, str,'location','northwest');
% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
%fig_pos = fig.PaperPosition;
%fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path1 = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = 'FS';
print(fig,[path1 fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[path1 fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path1 num2str(5)],'-dpng','-r600')
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file