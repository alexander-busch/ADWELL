% Model comparison


close all;
clear all;
clc;
addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange';
fig_path = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\LaTEX\Thesis\Figures\';

CreateFigure('Flowcurves', '','', {'lin' 'lin'}, '[1 1 21.0 7.4]');

subplot(1,2,1); hold on;
xlabel('Shear Rate [1/s]');
ylabel('Shear stress [Pa]');
set(gca,...
    'Xlim',[0 1200],...
    'Ylim',[0 16]);
box on;
grid on;

% Fann data
DR_Fann = [3 6 100 200 300 600];
RS_Fann = [6.57 6.85 12.20 15.01 18.29 28.14];

SR_Fann = 1.703*DR_Fann;
SS_Fann = 1.066*0.4788026*RS_Fann;

scatter(SR_Fann,SS_Fann);


% Bingham
PV = RS_Fann(6)-RS_Fann(5);
YP = RS_Fann(5)-PV;

SR = logspace(-2,4);
SS_Y = 1.066*0.4788026* (YP+PV/(1021.8-510.9)*SR);

plot(SR,SS_Y);

% Ostwald/PL
n = 3.32*log10(SS_Fann(6)/SS_Fann(5));
K = SS_Fann(5)/(511^n);

SS_PL = K.*SR.^n;

plot(SR,SS_PL);


% Herschel-Bulkley/YPL
tau_y = 2*SS_Fann(1)-SS_Fann(2);
n = 3.32*log10((SS_Fann(6) - tau_y)/(SS_Fann(5)-tau_y));
K = (SS_Fann(5)-tau_y)/(511^n);

SS_YPL = tau_y + K.*SR.^n;

plot(SR,SS_YPL);


% Cross model
eta_Fann = SS_Fann./SR_Fann;
[xData, yData] = prepareCurveData( SR_Fann, eta_Fann );
ft = fittype( 'mu_inf+(mu_0-mu_inf)/(1+(lambda*x)^(1-n))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf 0]; % opts.Lower = [0 4 0.001002 0];
opts.StartPoint = [0.63235924622541 0.0975404049994095 0.278498218867048 0.546881519204984];
opts.Upper = [1 +Inf +Inf 1]; % opts.Upper = [1 4 0.001002 1];
[fitresult, gof] = fit( xData, yData, ft, opts );

SS_Cross = (fitresult.mu_inf+(fitresult.mu_0-fitresult.mu_inf)./(1+(fitresult.lambda.*SR).^(1-fitresult.n))).*SR;

plot(SR,SS_Cross);

% (Caenn et al. 2011)
legend('Fann data',...
    'YP/PV (API 13D)',...
    'PL (API 13D)',...
    'YPL (API 13D)',...
    'Cross (1965)',...
    'location','southeast');

%TightFigure(2);


subplot(1,2,2); hold on;
xlabel('Shear Rate [1/s]');
ylabel('Apparent viscosity [Pa.s]');
set(gca,...
    'Xlim',[1e-2 10000],...
    'Ylim',[5e-3 1e2],...
    'XScale','log',...
    'YScale','log');
box on;
grid on;


% CreateFigure('Flowcurves', 'Shear Rate [1/s]','Apparent viscosity [Pa.s]', {'log' 'log'}, 'DINA5');


scatter(SR_Fann,SS_Fann./SR_Fann);
plot(SR,SS_Y./SR);
plot(SR,SS_PL./SR);
plot(SR,SS_YPL./SR);
plot(SR,SS_Cross./SR);

% legend('Fann data',...
%     'YP/PV (API 13D)',...
%     'PL (API 13D)',...
%     'YPL (API 13D)',...
%     'Cross (1965)');
%TightFigure(2);


% Print to files
fig_name = 'flowcurves';
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,fig_name,'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'





