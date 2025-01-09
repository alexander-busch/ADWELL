% Example

close all;
clear all;
clc;

% Fann readings
% Skalle (2013), page 26ff
% Drilling Fluids Processing Handbook, page 37ff
% Caenn et al. (2011), page 238ff
RS = [600 300	200 100	6 3]';
DR =  [140 98 78 54 16 13;...
    50 30 nan 20 8 nan;...
    30 19.5 16 13 7.3 7]';


%% Dial Reading vs. Rotor Speed
figure; hold on;
plot(RS,DR,'o');

% Bingham model based on Fann readings
PV = DR(1,:)-DR(2,:); % [cP = 1e-3 Pa.s]
YP = DR(2,:) - PV;  % [°]
tau_YPPV = YP+PV/300.*RS;
set(gca, 'ColorOrderIndex', 1);
plot(RS,tau_YPPV ,'-');

% PL model based on Fann readings
n_P = 3.32*log10((2*PV+YP)./(PV+YP));
n_P = 3.32*log10(DR(1,:)./DR(2,:));
K_P = (PV+YP)./(300.^n_P);
K_P = DR(2,:)./(300.^n_P);
tau_PL = K_P.*RS.^n_P;
set(gca, 'ColorOrderIndex', 1);
plot(RS,tau_PL,'--');


%% Dial Reading vs. Rotor Speed --> SI
figure; hold on;
plot(1.703*RS,0.4788026*1.066*DR,'o');
set(gca, 'ColorOrderIndex', 1);
plot(1.703*RS,0.4788026*1.066*tau_YPPV,'-');
set(gca, 'ColorOrderIndex', 1);
plot(1.703*RS,0.4788026*1.066*tau_PL,'--');


%% SI
% Bingham model based on SI values
SR = 1.703*RS;
PV = PV/1000; % [Pa.s]
YP = 0.4788026*1.066*YP;  % [Pa]
tau_YPPV = YP+PV.*(1.703*RS);
figure; hold on;
plot(SR,0.4788026*1.066*DR,'o');
set(gca, 'ColorOrderIndex', 1);
plot(SR,tau_YPPV ,'-');

% PL model based on SI values
n_P = 3.32*log10((0.4788026*1.066*DR(1,:))./(0.4788026*1.066*DR(2,:)));
K_P = 0.4788026*1.066*DR(2,:)./(511.^n_P);
tau_PL = K_P.*SR.^n_P;
set(gca, 'ColorOrderIndex', 1);
plot(SR,tau_PL,'--');




%% Sample calculation of pressure drop in drill pipe and annulus
% Caenn et al. (2011), page 238ff

close all;
clear all;
clc;
addpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic');

SR_max = 10000;
SR = logspace(0,log10(SR_max));
% Fann readings
rpm = [600 300	200	100	6 3]';
DR =  [30 19.5 16 13 7.3 7]';
% Conversion
SR_Fann = 1.703*rpm;
tau_Fann = 0.4788026*1.066*DR;
% Bingham model based on Fann readings
PV = DR(1)-DR(2); % [cP = 1e-3 Pa.s]
YP = DR(2) - PV;  % [lbf/100ft²]
tau_YPPV = 0.4788026*1.066*(YP+PV/(1022-511)*SR);

% Check dimensions of tau_0 and mu 
% tau_0 = 0.4788026*1.066*YP;
% mu = 0.4788026*1.066*PV/(1022-511);
% mu = PV/1000;
% tau = (tau_0 + mu*SR);
% plot(SR,tau);
% set(gca,'XScale','log','YScale','log');
% set(gca,'XScale','lin','YScale','lin');

% PL model based on six speed Fann readings (Pipe, Annulus)
n_P = 3.32*log10((2*PV+YP)/(PV+YP));
K_P = (PV+YP)/(511^n_P);
n_A6 = 0.657*log10(DR(4)/DR(6));
K_A6 = DR(4)/(170.3^n_A6);
tau_PL_P = 0.4788026*1.066*(K_P*SR.^n_P);
tau_PL6_A = 0.4788026*1.066*(K_A6*SR.^n_A6);
% PL model based on two speed Fann readings (Pipe, Annulus)
DR_100 = (YP+PV/(1022-511)*1.703*100);
DR_003 = (YP+PV/(1022-511)*1.703*003);
n_A2 = 0.657*log10(DR_100/DR_003);
K_A2 = DR_100/(170.3^n_A2);
tau_PL2_A = 0.4788026*1.066*(K_A2*SR.^n_A2);

fig1 = CreateFigure( 'Shear stress as a function shear rate',...
    'Shear rate [1/s]', 'Shear stress [Pa]', {'lin', 'lin'}, 'DINA5');
set(gca,...
    'Xlim',[0 1200],...
    'Ylim',[0 ceil(0.4788026*1.066*(YP+PV/(1022-511)*1200))]);
% Shaded shear rate intervalls for pipe and annular flow
SRI = area([SR_Fann(4) SR_Fann(6)],[SR_max SR_max]);
SRI.FaceColor = [0.8 0.8 0.8];
SRI.EdgeColor = 'none';
SRI.FaceAlpha = 0.5;
SRI = area([SR_Fann(1) SR_Fann(2)],[SR_max SR_max]);
SRI.FaceColor = [0.8 0.8 0.8];
SRI.EdgeColor = 'none';
SRI.FaceAlpha = 0.5;
% Plot data
plothandle=zeros(5,1);
plothandle(1) = plot(SR_Fann,tau_Fann,'o');
plothandle(2) = plot(SR,tau_YPPV);
plothandle(3) = plot(SR,tau_PL_P);
plothandle(4) = plot(SR,tau_PL6_A);
plothandle(5) = plot(SR,tau_PL2_A);
legendhandle = legend(plothandle,'Fann data',...
    'YP/PV (Bingham) model based on Fann dial readings at rpm = {300, 600}',...
    'PL (Ostwald) model based on Fann dial readings at rpm = {300, 600}',...
    'PL (Ostwald) model based on Fann dial readings at rpm = {3, 100}',...
    'PL (Ostwald) model based on YP/PV (Bingham) model at rpm = {3, 100}',...
    'Location','southeast');
% Plot text
% colors = get(groot,'defaultAxesColorOrder');
tx = 0.5*(SR_Fann(1)+SR_Fann(2));
ty = 15; % 0.4788026*1.066*(YP+PV/(1022-511)*tx)+0.1;
txt_P = text(tx,ty,...
    'Pipe flow',...
    'Color','w',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'Rotation',0);
tx = 0.5*(SR_Fann(4)+SR_Fann(6));
% ty = 0.4788026*1.066*(YP+PV/(1022-511)*tx)+0.5;
txt_A = text(tx,ty,...
    {'Annular','flow'},...
    'Color','w',...
    'FontSize',14,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'Rotation',0);
% Second set of axis
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Xlim',xlim/1.703,...
    'Ylim',ylim/0.4788026,...
    'Color','none');
xlabel('Fann viscometer speed [rpm]');
ylabel('Shear stress [lbf/100ft²]');

% Copy figure --> Shear stress (lin)
fig2 = figure;
set(fig2,...
    'Name','Shear stress as a function shear rate',...
    'color','w',...
    'Units','centimeters',...
    'Position',[1 1 21.0 14.8]); % DIN A5 ;
objects=allchild(fig1);
copyobj(get(fig1,'children'),fig2);

% Modify original figure --> (log)
set([ax1 ax2],'XScale','log','YScale','log');
set(txt_P,...
    'String',{'Pipe','flow'},...
    'Position',[730 2]);
set(txt_A,...
    'String','Annular flow',...
    'Position',[30 10]);

% Copy figure --> Shear stress (log)
fig3 = figure;
set(fig3,...
    'Name','Shear stress as a function shear rate',...
    'color','w',...
    'Units','centimeters',...
    'Position',[1 1 21.0 14.8]); % DIN A5 ;
copyobj(get(fig1,'children'),fig3);

% Modify original figure --> Apparent viscosity (log)
fig1;
set(fig1,'Name','Apparent viscosity as a function shear rate');
set(plothandle(1),'YData',tau_Fann./SR_Fann);
set(plothandle(2),'YData',tau_YPPV./SR);
set(plothandle(3),'YData',tau_PL_P./SR);
set(plothandle(4),'YData',tau_PL6_A./SR);
set(plothandle(5),'YData',tau_PL2_A./SR);
set(legendhandle,'Location','southwest');
ylabel(ax1,'Apparent viscosity [Pa.s]');
ylabel(ax2,'Apparent viscosity [cP]');
set(ax1,'Ylim',[1e-3 1e1]);
set(ax2,'Ylim',[1e-3 1e1]/10);
set(txt_A,'Position',[30 2]);

% Wall shear rate
SR_W_P = (1+3*n_P)/(4*n_P)*SR;
SR_W_A6 = (1+2*n_A6)/(3*n_A6)*SR;
SR_W_A2 = (1+2*n_A2)/(3*n_A2)*SR;
% figure;
% hold on;
% plot(SR,SR_W_P./SR);
% plot(SR,SR_W_A6./SR);
% plot(SR,SR_W_A2./SR);