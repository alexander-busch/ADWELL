clear all;
close all;
clc;

col_PAC8 = [0,0,1]; 
col_PAC4 = [0,0.5,1]; 
col_PAC2 = [0,1,1];
special = [0,0,0,]; %[34, 177, 76]/255;

addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\stress-based';

load FC
load FC_Cross_PAC2;
load FC_Cross_PAC4;
load FS_PAC2;
load FS_PAC4;

% Define data range for PL-fit to FS data
DataRange_PAC2 = [1:length(AngularFrequency_PAC2)];
DataRange_PAC4 = [1:length(AngularFrequency_PAC4)];


%% Figure

fig1 = figure('color','w','Units','centimeters','Position',[1 1 21 14.8]); % DIN A5 
hold on;
grid on;
box on;

yyaxis left;
xlabel('Shear rate [1/s], Angular Frequency [rad/s]');
ylabel('Shear stress [Pa], FNSD [Pa]');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'YColor','k',...
    'xlim', [1e0 1e2],...
    'ylim', [1e-2 2e2]);

yyaxis right;
ylabel('Recoverable shear [-]');
set(gca,...
    'YScale','log',...
    'YColor','k',...
    'ylim', [1e-1 1e1]);


%% PAC2
yyaxis left;
% Shear stress
ApparentViscosity=mu_inf_PAC2+(mu_0_PAC2-mu_inf_PAC2)./(1+(lambda_Cr_PAC2.*ShearRate_PAC2).^n_Cr_PAC2);
ShearStress = ApparentViscosity.*ShearRate_PAC2;
g(1) = plot(ShearRate_PAC2,ShearStress,'-','Color',col_PAC2);

% First normal stress difference estimated with Launs rule
FirstNormalStressCoeff_PAC2 = 2.*(StorageModulus_PAC2./(AngularFrequency_PAC2.^2)).*(1+(StorageModulus_PAC2./LossModulus_PAC2).^2).^0.5;
FirstNormalStressDiff_PAC2 = FirstNormalStressCoeff_PAC2.*AngularFrequency_PAC2.^2;
g(2) = plot(AngularFrequency_PAC2,FirstNormalStressDiff_PAC2,'d','Color',col_PAC2,'MarkerFaceColor',col_PAC2);

% Piecewise PL-fit of FNSD
PL_FNSD_PAC2_high = createFit_PowerLaw(AngularFrequency_PAC2(1:4),FirstNormalStressDiff_PAC2(1:4));
close gcf;
PL_FNSD_PAC2_low = createFit_PowerLaw(AngularFrequency_PAC2(5:end),FirstNormalStressDiff_PAC2(5:end));
close gcf;

ShearRate_PAC2_high = linspace( 25.1,100);
ShearRate_PAC2_low = linspace(1, 25.1);

FNSD_fit_high=PL_FNSD_PAC2_high.K.*ShearRate_PAC2_high.^(PL_FNSD_PAC2_high.n);
FNSD_fit_low=PL_FNSD_PAC2_low.K.*ShearRate_PAC2_low.^(PL_FNSD_PAC2_low.n);
g(3) = plot(ShearRate_PAC2_high,FNSD_fit_high,'--','Color',col_PAC2,'LineWidth',1);
g(3) = plot(ShearRate_PAC2_low,FNSD_fit_low,'--','Color',col_PAC2,'LineWidth',1);


% Stress ratio & recoverable shear
yyaxis right;
ApparentViscosity=mu_inf_PAC2+(mu_0_PAC2-mu_inf_PAC2)./(1+(lambda_Cr_PAC2.*ShearRate_PAC2_high).^n_Cr_PAC2);
ShearStress = ApparentViscosity.*ShearRate_PAC2_high;
StressRatio = ( PL_FNSD_PAC2_high.K.*ShearRate_PAC2_high.^(PL_FNSD_PAC2_high.n) )./ShearStress;
% g(4) = plot(ShearRate_PAC2,StressRatio,'k-','LineWidth',2);
g(4) = plot(ShearRate_PAC2_high,StressRatio/2,'-','Color',col_PAC2,'LineWidth',2);

ApparentViscosity=mu_inf_PAC2+(mu_0_PAC2-mu_inf_PAC2)./(1+(lambda_Cr_PAC2.*ShearRate_PAC2_low).^n_Cr_PAC2);
ShearStress = ApparentViscosity.*ShearRate_PAC2_low;
StressRatio = ( PL_FNSD_PAC2_low.K.*ShearRate_PAC2_low.^(PL_FNSD_PAC2_low.n) )./ShearStress;
% g(4) = plot(ShearRate_PAC2,StressRatio,'k-','LineWidth',2);
g(4) = plot(ShearRate_PAC2_low,StressRatio/2,'-','Color',col_PAC2,'LineWidth',2);




%% PAC4
yyaxis left;
% Shear stress
ApparentViscosity=mu_inf_PAC4+(mu_0_PAC4-mu_inf_PAC4)./(1+(lambda_Cr_PAC4.*ShearRate_PAC4).^n_Cr_PAC4);
ShearStress = ApparentViscosity.*ShearRate_PAC4;
g(5) = plot(ShearRate_PAC4,ShearStress,'-','Color',col_PAC4);

% First normal stress difference estimated with Launs rule
FirstNormalStressCoeff_PAC4 = 2.*(StorageModulus_PAC4./(AngularFrequency_PAC4.^2)).*(1+(StorageModulus_PAC4./LossModulus_PAC4).^2).^0.5;
FirstNormalStressDiff_PAC4 = FirstNormalStressCoeff_PAC4.*AngularFrequency_PAC4.^2;
g(6) = plot(AngularFrequency_PAC4,FirstNormalStressDiff_PAC4,'d','Color',col_PAC4,'MarkerFaceColor',col_PAC4);

% Piecewise PL-fit of FNSD
PL_FNSD_PAC4_high = createFit_PowerLaw(AngularFrequency_PAC4(1:4),FirstNormalStressDiff_PAC4(1:4));
close gcf;
PL_FNSD_PAC4_low = createFit_PowerLaw(AngularFrequency_PAC4(5:end),FirstNormalStressDiff_PAC4(5:end));
close gcf;

ShearRate_PAC4_high = linspace( 25.1,100);
ShearRate_PAC4_low = linspace(1, 25.1);

FNSD_fit_high=PL_FNSD_PAC4_high.K.*ShearRate_PAC4_high.^(PL_FNSD_PAC4_high.n);
FNSD_fit_low=PL_FNSD_PAC4_low.K.*ShearRate_PAC4_low.^(PL_FNSD_PAC4_low.n);
g(7) = plot(ShearRate_PAC4_high,FNSD_fit_high,'--','Color',col_PAC4,'LineWidth',1);
g(7) = plot(ShearRate_PAC4_low,FNSD_fit_low,'--','Color',col_PAC4,'LineWidth',1);


% Stress ratio & recoverable shear
yyaxis right;
ApparentViscosity=mu_inf_PAC4+(mu_0_PAC4-mu_inf_PAC4)./(1+(lambda_Cr_PAC4.*ShearRate_PAC4_high).^n_Cr_PAC4);
ShearStress = ApparentViscosity.*ShearRate_PAC4_high;
StressRatio = ( PL_FNSD_PAC4_high.K.*ShearRate_PAC4_high.^(PL_FNSD_PAC4_high.n) )./ShearStress;
% g(4) = plot(ShearRate_PAC4,StressRatio,'k-','LineWidth',2);
g(8) = plot(ShearRate_PAC4_high,StressRatio/2,'-','Color',col_PAC4,'LineWidth',2);

ApparentViscosity=mu_inf_PAC4+(mu_0_PAC4-mu_inf_PAC4)./(1+(lambda_Cr_PAC4.*ShearRate_PAC4_low).^n_Cr_PAC4);
ShearStress = ApparentViscosity.*ShearRate_PAC4_low;
StressRatio = ( PL_FNSD_PAC4_low.K.*ShearRate_PAC4_low.^(PL_FNSD_PAC4_low.n) )./ShearStress;
% g(4) = plot(ShearRate_PAC4,StressRatio,'k-','LineWidth',2);
g(8) = plot(ShearRate_PAC4_low,StressRatio/2,'-','Color',col_PAC4,'LineWidth',2);




%% Legend
h_legend = legend( g(1:8),...
    'PAC2 - Shear stress $\tau (\dot \gamma)$, Cross fit of FC data', ...
    'PAC2 - FNSD $N_{1} = \psi_{1} \dot \gamma^2$, $\psi_{1}$ obtained from FS data (Laun 1986)', ...
    ['PAC2 - FNSD $N_{1}$, PL fit (m'' = ' num2str(PL_FNSD_PAC2.K,3) ', n'' = ' num2str(PL_FNSD_PAC2.n,3) ')' ], ...
    'PAC2 - Recoverable shear $\frac{N_{1}}{2 \tau}$[-]',...
    'PAC4 - Shear stress $\tau (\dot \gamma)$, Cross fit of FC data', ...
    'PAC4 - FNSD $N_{1} = \psi_{1} \dot \gamma^2$, $\psi_{1}$ obtained from FS data (Laun 1986)', ...
    ['PAC4 - FNSD $N_{1}$, PL fit (m'' = ' num2str(PL_FNSD_PAC4.K,3) ', n'' = ' num2str(PL_FNSD_PAC4.n,3) ')' ], ...
    'PAC4 - Recoverable shear $\frac{N_{1}}{2 \tau}$[-]');

set(h_legend,...
    'Location','northwest',...
    'Interpreter','Latex');


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
set(fig1,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig1.PaperPosition;
fig1.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub1\Figures\';
fig_name = 'StressRatio';
print(fig1,[path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig1,[path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig1,[path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file