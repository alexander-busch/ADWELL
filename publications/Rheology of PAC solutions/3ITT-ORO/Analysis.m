%% Clean up
clear all;
close all;
clc;

%% Import SP data

[~, ~, raw] = xlsread([pwd '\3ITT-ORO.xlsx'],'SP','A3:M160');

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));
% Allocate imported array to column variable names
Time = data(:,1);
AngularFrequency = data(:,2);
Strain = data(:,3);
ShearRate = data(:,4);
ShearStress = data(:,5);
Viscosity = data(:,6);
ComplexModulus = data(:,7);
StorageModulus = data(:,8);
LossModulus = data(:,9);
DampingFactor = data(:,10);
ComplexViscosity = data(:,11);
DeflectionAngle = data(:,12);
Torque = data(:,13);

% Clear temporary variables
clearvars data raw R;
%% Create figure
fig = figure('color','w','Units','centimeters','Position',[1 1 21 14.8]); % DIN A5 

MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor
hold on;


%% Formating
ylim_min = 0;
ylim_max = 1.2;
xlabel('Time [s]');
ylabel('G'', G'''' [Pa]');
grid('on');
set(gca,...
    'xlim', [0 202],... % Max of SP data = 1530
    'ylim', [ylim_min ylim_max],... 
    'box','on');



% Interval 2
HSR = area([100 101], [ylim_max ylim_max]);
HSR.FaceColor = [0.8 0.8 0.8];
HSR.EdgeColor = 'none';
HSR.FaceAlpha = 0.5;


% Interval shear rates labeling
y1 = 0.6;
dy1= 0.1;
dy2=0.06;

x1 = 50;
txt1 = '($\omega~=~10~rad/s,~\gamma~=~0.1~\%$)';
text(x1,y1,txt1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'interpreter','latex');
text(x1,y1+dy2,'Oscillation','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);

x1 = 100.5;
txt1 = '($\dot{\gamma}~=~3000~s^{-1}$)';
text(x1,y1+dy1,txt1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'interpreter','latex');
text(x1,y1+dy1+dy2,'Rotation','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);

x1 = 150;
txt1 = '($\omega~=~10~rad/s,~\gamma~=~0.1~\%$)';
text(x1,y1,txt1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'interpreter','latex');
text(x1,y1+dy2,'Oscillation','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);

%% Plot SP data

MFC = 'k'; % MarkerFaceColor

g(1) = plot(Time,LossModulus,...
    'kd','markersize',MS,'MarkerFaceColor',MFC);
g(2) = plot(Time,StorageModulus,...
    'ks','markersize',MS,'MarkerFaceColor',MFC);


%% Additional information

plot([0 100], [9.53e-1 9.53e-1],'k-');
plot([101 1400], [9.39e-1 9.39e-1],'k-');
plot([0 100], [4.57e-1 4.57e-1],'k-');
plot([101 1400], [2.37e-1 2.37e-1],'k-');

plot([0 100], [1.01 1.01],'k--');
plot([101 1400], [9.943e-1 9.943e-1],'k--');
plot([0 100], [2.522e-1 2.522e-1],'k--');
plot([101 1400], [2.797e-1 2.797e-1],'k--');




%% Import UiS data
[~, ~, raw] = xlsread([pwd '\3ITT-ORO.xlsx'],'UiS','A3:M113');

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));
% Allocate imported array to column variable names
Time = data(:,1);
AngularFrequency = data(:,2);
Strain = data(:,3);
ShearRate = data(:,4);
ShearStress = data(:,5);
Viscosity = data(:,6);
ComplexModulus = data(:,7);
StorageModulus = data(:,8);
LossModulus = data(:,9);
DampingFactor = data(:,10);
ComplexViscosity = data(:,11);
DeflectionAngle = data(:,12);
Torque = data(:,13);

% Clear temporary variables
clearvars data raw R;


%% Plot UiS data

MFC = 'w'; % MarkerFaceColor

g(3) = plot(Time,LossModulus,...
    'kd','markersize',MS,'MarkerFaceColor',MFC);
g(4) = plot(Time,StorageModulus,...
    'ks','markersize',MS,'MarkerFaceColor',MFC);


%% Legend

legend( g(1:4),...
        'Lab1 - Loss modulus G''''',...
        'Lab1 - Storage modulus G''',...
        'Lab2 - Loss modulus G''''',...
        'Lab2 - Storage modulus G''');
%         'Location','southeast',...
%         'FontSize',21);




% text(25,0.43,'Storage modulus G''','HorizontalAlignment','center','FontSize',10);
% text(25,0.93,'Loss modulus G''''','HorizontalAlignment','center','FontSize',10);


% h_legend = legend(...
%     'Loss modulus G''''', ...
%     'Storage modulus G''');
% 
% set(h_legend,'Location','southwest','FontSize',21);


%% Print as image

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

% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Figures\';
fig_name = '3ITT-ORO';
print(fig,[path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg'
print(fig,[path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_nam],'pdf') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file