% Analysis of distribution of app viscositites for every shear rate of Lab1
% Flowcurves

%% Create figure

close all;
fig = figure; % ('units','normalized','outerposition',[0 0 1 1]);
hold on;

%% Formating

xlabel(cat(2,headers{1,3},' ',headers{2,3}));
ylabel(cat(2,headers{1,4},' ',headers{2,4}));
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-2 1.2e3],...
    'ylim', [2e-2 3e-1],...
    'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1],...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');


%% Plot raw data of all Bergen flow curves (Lab1)

hold on;
for i = 1:13
    plot(PAC4flowcurves(:,3,i),PAC4flowcurves(:,4,i),'ko');
end

%% Plot raw data of all UiS flow curves (Lab2)

% Scatter of Milads data
load Lab2_Scatter2016;
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)+ 3*Lab2_Scatter2016(:,3),'r--');
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)- 3*Lab2_Scatter2016(:,3),'r--');


%% Create PDF - Upward sweep

% shearrate = PAC4flowcurves(40,3,1);
% 0.1 < gamma_dot < 100 - 21:40
% 100 < gamma_dot < 1200 - 41:50

SR_idx_min = 1;
SR_idx_max = 50;


h_sum = zeros(1,24);
h_all = zeros(SR_idx_max - SR_idx_min-1,24);
fig = figure();
hold on;

for i = SR_idx_min:SR_idx_max % Upward 21:50 for shear rate > 0.1
    viscosity = PAC4flowcurves(i,4,:); % Get all viscosity values for respective shearateindex
    vis_Mean = sum(viscosity)/length(viscosity); % Determine mean of viscosity MATLAB mean produces error...

    viscosity = viscosity ./ vis_Mean - 1; % Non-dimensionalise viscosity values around zero
    vis_SD = std(viscosity); % Compute SD of non-dim. viscosities

    m3SD = -3*vis_SD;
    m2SD = -2*vis_SD;
    m1SD = -1*vis_SD;
    mean = 0;
    p1SD = 1*vis_SD;
    p2SD = 2*vis_SD;
    p3SD = 3*vis_SD;

    X_PDF = linspace(m3SD,p3SD,25); % Define non-dimensional SD vector from -3SD to + 3SD
    
    h = histogram(viscosity,X_PDF,'Visible','off');
    h_sum = h_sum + h.Values;
    h_all(i+1-SR_idx_min,:) = h.Values;
end

X_PDF = linspace(-3,3,24);
bar(X_PDF,h_sum,'FaceColor','k');

figure
b = bar3(X_PDF,h_all');
h=findobj('type','surface');
set(h,'facecolor',[0.7 0.7 0.7])

%Format 2D

xlabel('Non-dim. mean +/- 3SD');
ylabel('Cummulated frequencies [-]');

grid on
set(gca,...
    'xlim', [-3 3],...
    'XTick', [-3 -2 -1 0 1 2 3],...
    'box','on',...
    'XColor', 'k',...
    'YColor', 'k',...
    'FontSize',24);
set(gcf,...
    'color','w'); 

%Format 3D

xlabel('Shear rate index');
ylabel('Non-dim. mean +/- 3SD');
zlabel('Cummulated frequencies [-]');

grid on
set(gca,...
    'ylim', [-3 3],...
    'yTick', [-3 -2 -1 0 1 2 3],...
    'box','on',...
    'XColor', 'k',...
    'YColor', 'k',...
    'ZColor', 'k',...
    'FontSize',24);
set(gcf,...
    'color','w'); 



%% Create PDF - Downward sweep

% shearrate = PAC4flowcurves(91,3,1);
% 100 > gamma_dot > 0.1 - 72:90
% 1200 > gamma_dot > 100 - 61:71

SR_idx_min = 61;
SR_idx_max = 110;

h_sum = zeros(1,24);
h_all = zeros(SR_idx_max - SR_idx_min-1,24);
fig = figure();
hold on;

for i = SR_idx_min:SR_idx_max % Downward 61:80 for shear rate > 0.1
    viscosity = PAC4flowcurves(i,4,:); % Get all viscosity values for respective shearateindex
    vis_Mean = sum(viscosity)/length(viscosity); % Determine mean of viscosity MATLAB mean produces error...

    viscosity = viscosity ./ vis_Mean - 1; % Non-dimensionalise viscosity values around zero
    vis_SD = std(viscosity); % Compute SD of non-dim. viscosities

    m3SD = -3*vis_SD;
    m2SD = -2*vis_SD;
    m1SD = -1*vis_SD;
    mean = 0;
    p1SD = 1*vis_SD;
    p2SD = 2*vis_SD;
    p3SD = 3*vis_SD;

    X_PDF = linspace(m3SD,p3SD,25); % Define non-dimensional SD vector from -3SD to + 3SD
    
    h = histogram(viscosity,X_PDF,'Visible','off');
    % bar3(x,h.Values);
    h_sum = h_sum + h.Values;
end

X_PDF = linspace(-3,3,24);
bar(X_PDF,h_sum,'FaceColor','k');

