clear;
clc;
close all;

path1 = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\MATLAB\AdWell\Rheology\PAC\ParameterSpace\';
path2 = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\MATLAB\AdWell\Rheology\PAC\ParameterSpace\';

%% Import data from spreadsheet

% Import the AS data
[~, ~, raw] = xlsread( [path1 'AS_ParameterRange.xlsx'],'Sheet1','A2:C27');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
AS_omega = data(:,1);
AS_gamma = data(:,2);
AS_gamma_dot = data(:,3);

% Clear temporary variables
clearvars data raw;

% Import the FS data
[~, ~, raw] = xlsread( [path2 'FS_ParameterRange.xlsx'],'Sheet1','A2:C20');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
FS_omega = data(:,1);
FS_gamma = data(:,2);
FS_gamma_dot = data(:,3);

% Clear temporary variables
clearvars data raw;

%% 


figure;
hold on;
plot(AS_gamma, AS_omega,'ko');
plot(FS_gamma, FS_omega, 'ko');


% Lines of constant shear rate
axis = logspace(-3,4,8);
plot(axis, flip(axis)); % Axis

for i = 1:length(axis)
    for j = 1:length(axis)
        text(axis(i),axis(j),num2str(axis(i)/100*axis(j)));
        plot(axis(i), flip(axis(j)));
    end;
end;

AS_gamma_dot(1)
AS_gamma_dot(length(AS_gamma_dot))









xlabel('Strain Gamma [%]');
ylabel('Angular frequency Omega [rad/s]');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [min(axis) max(axis)],...
    'ylim', [min(axis) max(axis)],...
    'XTick', axis,...
    'YTick', axis,...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');
