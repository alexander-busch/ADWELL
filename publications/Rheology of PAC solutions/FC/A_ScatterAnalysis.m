%% Analysis of PAC4 scatter study data
% Rheometric data obtained at SINTEF Petroleaum AS, Bergen, 09.07.15
% (C) by Alexander Busch, NTNU, Dec. 2015

clc;
clear;
close all;


% Data
path_lab1 = [pwd '\data\Lab1_SINTEF_Petroleum\'];
path_lab2 = [pwd '\data\Lab2_UiS\'];
path_lab3 = [pwd '\data\Lab3_NTNU\'];
path_Literature = [pwd '\data\Literature\'];

addpath ([pwd '\functions\']);
addpath (path_lab1);
addpath (path_lab2);
addpath (path_lab3);
addpath (path_Literature);

%% Import of data

PAC2data = importdata([path_lab1 'PAC2.csv']);
PAC4data = importdata([path_lab1 'PAC4.csv']);
PAC8data = importdata([path_lab1 'PAC8.csv']);

headers = importheaders([path_lab1 'PAC4.csv']);

% Raw data intervall specification
Up = [1,50]; % Upward flow curve shear rate range
Max = [62,71]; % Shear rate = 1200 1/s
Down = [83,132]; % Downward flow curve shear rate range
intervall = 165; % Intervall of data series in CSV

% Determination of imported data intervalls
Up = [1,50];
Max = [Up(2)+1, Up(2)+1+(Max(2)-Max(1))];
Down = [Max(2)+1, Max(2)+1+(Down(2)-Down(1))];   


%% Creation of 3D data matrices

PAC2flowcurves = cat(1, PAC2data(1:50,1:4), PAC2data(62:71,1:4), PAC2data(83:132,1:4)); % First Flowcurve matrix
    sp = intervall;
    flowcurve = cat(1, PAC2data(sp:sp+49,1:4), PAC2data(sp+61:sp+70,1:4), PAC2data(sp+82:sp+131,1:4));
    PAC2flowcurves = cat(3, PAC2flowcurves, flowcurve);

PAC4flowcurves = cat(1, PAC4data(1:50,1:4), PAC4data(62:71,1:4), PAC4data(83:132,1:4)); % First Flowcurve matrix
for i = 1:12 % Creation of second to thirteenth Flowcurve matrix and final 3D matrix
    sp = i*intervall-(i-1);
    flowcurve = cat(1, PAC4data(sp:sp+49,1:4), PAC4data(sp+61:sp+70,1:4), PAC4data(sp+82:sp+131,1:4));
    PAC4flowcurves = cat(3, PAC4flowcurves, flowcurve);
end

PAC8flowcurves = cat(1, PAC8data(1:50,1:4), PAC8data(62:71,1:4), PAC8data(83:132,1:4)); % First Flowcurve matrix
    sp = intervall;
    flowcurve = cat(1, PAC8data(sp:sp+49,1:4), PAC8data(sp+61:sp+70,1:4), PAC8data(sp+82:sp+131,1:4));
    PAC8flowcurves = cat(3, PAC8flowcurves, flowcurve);


%% Statistical analysis

PAC2_SS_SD = std(PAC2flowcurves(:,2,:),0,3); % Shear stress standard deviation
PAC2_SS_M  = mean(PAC2flowcurves(:,2,:),3); % Shear stress mean
PAC2_AV_SD = std(PAC2flowcurves(:,4,:),0,3); % Apparent viscosity standard deviation
PAC2_AV_M  = mean(PAC2flowcurves(:,4,:),3); % Apparent viscosity standard mean

PAC4_SS_SD = std(PAC4flowcurves(:,2,:),0,3); % Shear stress standard deviation
PAC4_SS_M  = mean(PAC4flowcurves(:,2,:),3); % Shear stress mean
PAC4_AV_SD = std(PAC4flowcurves(:,4,:),0,3); % Apparent viscosity standard deviation
PAC4_AV_M  = mean(PAC4flowcurves(:,4,:),3); % Apparent viscosity standard mean

PAC8_SS_SD = std(PAC8flowcurves(:,2,:),0,3); % Shear stress standard deviation
PAC8_SS_M  = mean(PAC8flowcurves(:,2,:),3); % Shear stress mean
PAC8_AV_SD = std(PAC8flowcurves(:,4,:),0,3); % Apparent viscosity standard deviation
PAC8_AV_M  = mean(PAC8flowcurves(:,4,:),3); % Apparent viscosity standard mean

%% Prepare vectors to plot

% Shear rate
Up_SR = PAC4flowcurves(Up(1):Up(2),3,1); % Up
Down_SR = PAC4flowcurves(Down(1):Down(2),3,1); % Down

% Shear stress
Up_SS_M = PAC4_SS_M(Up(1):Up(2));
Up_SS_M_plus3D = PAC4_SS_M(Up(1):Up(2))+3.*PAC4_SS_SD(Up(1):Up(2));
Up_SS_M_minus3D = PAC4_SS_M(Up(1):Up(2))-3.*PAC4_SS_SD(Up(1):Up(2));
Down_SS_M = PAC4_SS_M(Down(1):Down(2));
Down_SS_M_plus3D = PAC4_SS_M(Down(1):Down(2))+3.*PAC4_SS_SD(Down(1):Down(2));
Down_SS_M_minus3D = PAC4_SS_M(Down(1):Down(2))+3.*PAC4_SS_SD(Down(1):Down(2));

% Apparent viscosity
Up_AV_M = PAC4_AV_M(Up(1):Up(2));
Up_AV_M_plus3D = PAC4_AV_M(Up(1):Up(2))+3.*PAC4_AV_SD(Up(1):Up(2));
Up_AV_M_minus3D = PAC4_AV_M(Up(1):Up(2))-3.*PAC4_AV_SD(Up(1):Up(2));
Down_AV_M = PAC4_AV_M(Down(1):Down(2));
Down_AV_M_plus3D = PAC4_AV_M(Down(1):Down(2))+3.*PAC4_AV_SD(Down(1):Down(2));
Down_AV_M_minus3D = PAC4_AV_M(Down(1):Down(2))-3.*PAC4_AV_SD(Down(1):Down(2));

parentpath = cd(cd('..'));
path = [parentpath '\Comparison_FC_3ITT'];
save([path '\FC_Lab1.mat'],...
    'Up_SR',...
    'Down_SR',...
    'Up_AV_M',...
    'Down_AV_M');




%% Plot

% Axes
AxesDefinition = 'log'; % normal or log

% Stress vs. Shear Rate (All flow curves)
figure('Name','Stress vs. Shear Rate (All flow curves)','color','w','NumberTitle','off');
grid('on');
hold on;
for i = 1:13
    plot(PAC4flowcurves(:,3,i),PAC4flowcurves(:,2,i));
end
xlabel(cat(2,headers{1,3},' ',headers{2,3})), ylabel(cat(2,headers{1,2},' ',headers{2,2}));

% Viscosity vs. Shear Rate (All flow curves)
figure('Name','Viscosity vs. Shear Rate (All flow curves)','color','w','NumberTitle','off');
grid('on');
hold on;
for i = 1:13
    plot(PAC4flowcurves(:,3,i),PAC4flowcurves(:,4,i));
end
set(gca,'XScale','log','YScale','log');
xlabel(cat(2,headers{1,3},' ',headers{2,3})), ylabel(cat(2,headers{1,4},' ',headers{2,4}));

% Stress vs. Shear Rate (Mean & SD)
Stat_SR_SS = figure(...
    'Name','Stress vs. Shear Rate (Mean & SD)',...
    'color','w','NumberTitle','off');
hold on;
plot(Up_SR, Up_SS_M,'g','LineWidth',1);
plot(Up_SR, Up_SS_M_plus3D,'g--');
plot(Up_SR, Up_SS_M_minus3D,'g--');
plot(Down_SR, Down_SS_M,'LineWidth',1),'b';
plot(Down_SR, Down_SS_M_plus3D,'b--');
plot(Down_SR, Down_SS_M_minus3D,'b--');
title('');
xlabel(cat(2,headers{1,3},' ',headers{2,3})), ylabel(cat(2,headers{1,2},' ',headers{2,2}));
legend('Upward Mean', 'Upward Mean + 3SD', 'Upward Mean - 3SD', 'Downward Mean', 'Downward Mean + 3SD', 'Downward Mean - 3SD', 'Location','southeast');
grid('on');
set(gca,'XScale','log','YScale','log');

% Viscosity vs. Shear Rate (Mean & SD)
Stat_SR_AV = figure(...
    'Name','Viscosity vs. Shear Rate (Mean & SD)',...
    'color','w','NumberTitle','off');
hold on;
plot(Up_SR, Up_AV_M,'g','LineWidth',1);
plot(Up_SR, Up_AV_M_plus3D,'g--');
plot(Up_SR, Up_AV_M_minus3D,'g--');
plot(Down_SR, Down_AV_M,'LineWidth',1),'b';
plot(Down_SR, Down_AV_M_plus3D,'b--');
plot(Down_SR, Down_AV_M_minus3D,'b--');title('');

load Fann.mat; % Lab1 - Fann viscometer
plot (Fann(1,:),Fann(2,:)./Fann(1,:)+ 3*Fann(3,:)./Fann(1,:),'r--');
plot (Fann(1,:),Fann(2,:)./Fann(1,:)- 3*Fann(3,:)./Fann(1,:),'r--');

load Lab2_Scatter2016.mat; % Lab2 - UiS data
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)+ 3*Lab2_Scatter2016(:,3),'k+');
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)- 3*Lab2_Scatter2016(:,3),'k+');

xlabel(cat(2,headers{1,3},' ',headers{2,3})), ylabel(cat(2,headers{1,4},' ',headers{2,4}));
legend('Upward Mean', 'Upward Mean + 3SD', 'Upward Mean - 3SD', 'Downward Mean', 'Downward Mean + 3SD', 'Downward Mean - 3SD', 'Fann Mean + 3SD', 'Fann Mean - 3SD', 'UiS Mean + 3SD', 'UiS Mean - 3SD', 'Location','northeast');
grid('on');
if strcmp(AxesDefinition,'log') == 1
    set(gca,...
        'XScale','log','YScale','log',...
        'xlim', [1e-2 2e3],...
        'ylim', [2e-2 3e-1]);
end