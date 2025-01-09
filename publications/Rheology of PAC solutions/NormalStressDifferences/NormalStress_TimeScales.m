clear all;
close all;
clc;

col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 



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

fig2 = figure('color','w','Units','centimeters','Position',[1 1 12 8.45]); % DIN A5 
hold on;
grid on;
box on;

yyaxis left;
xlabel('$\dot\gamma$ [1/s], $\omega$ [rad/s]','FontSize',10,'Interpreter','latex');
ylabel('$\eta(\dot\gamma)$, $\eta^{*}(\omega=\dot\gamma)$ [Pa.s], $\psi_{1}$ [Pa.s$^2$]','FontSize',10,'Interpreter','latex');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'YColor','k',...
    'xlim', [1e0 1e2],...
    'ylim', [2e-3 2e-1]);

yyaxis right;
ylabel('$\lambda_{rs}$ [-]','FontSize',10,'Interpreter','latex');
set(gca,...
    'YScale','log',...
    'YColor','k',...
    'ylim', [2e-2 2e0]);


%% PAC2
yyaxis left;

% Apparent viscosity estimated, Cox-Merz rule
g(1) = plot(AngularFrequency_PAC2,ComplexViscosity_PAC2,'o','Color',col_PAC2,'MarkerFaceColor',col_PAC2);

% Apparent viscosity, Cross fit
ShearRate = [logspace(-1,3) 1200];
ApparentViscosity=mu_inf_PAC2+(mu_0_PAC2-mu_inf_PAC2)./(1+(lambda_Cr_PAC2.*ShearRate).^n_Cr_PAC2);
g(2) = plot(ShearRate,ApparentViscosity,'-','Color',col_PAC2);

% First normal stress coefficient, Launs rule
FirstNormalStressCoeff_PAC2 = 2.*(StorageModulus_PAC2./(AngularFrequency_PAC2.^2)).*(1+(StorageModulus_PAC2./LossModulus_PAC2).^2).^0.5;
FirstNormalStressDiff_PAC2 = FirstNormalStressCoeff_PAC2.*AngularFrequency_PAC2.^2;
%g(3) = plot(AngularFrequency_PAC2,FirstNormalStressCoeff_PAC2,'d','Color',col_PAC2,'MarkerFaceColor',col_PAC2);

% PL-fit of FNSD & FNSC
PL_FNSD_PAC2 = 	createFit_PowerLaw(AngularFrequency_PAC2(DataRange_PAC2),FirstNormalStressDiff_PAC2(DataRange_PAC2));
close gcf;
FNSC_fit = PL_FNSD_PAC2.K.*ShearRate.^(PL_FNSD_PAC2.n-2);
% Element-wise n and K with fit function of n and K
% [ FNSC_fit, K_coeff, n_coeff ] = get_FNSC_fit_variable_N_K( AngularFrequency_PAC2(DataRange_PAC2),FirstNormalStressDiff_PAC2(DataRange_PAC2),ShearRate );
g(3) = plot(ShearRate,FNSC_fit,'--','Color',col_PAC2);

% Compute local Power Law n & K for Cross fit of apparent viscosity data
k_PL = exp(log(ShearRate).*(mu_0_PAC2-mu_inf_PAC2).*n_Cr_PAC2.*lambda_Cr_PAC2.*ShearRate./((mu_inf_PAC2.*(lambda_Cr_PAC2.*ShearRate+1).^n_Cr_PAC2+mu_0_PAC2-mu_inf_PAC2).*(lambda_Cr_PAC2.*ShearRate+1))).*(mu_inf_PAC2+(lambda_Cr_PAC2.*ShearRate+1).^(-n_Cr_PAC2)*mu_0_PAC2-(lambda_Cr_PAC2.*ShearRate+1).^(-n_Cr_PAC2).*mu_inf_PAC2);
n_PL=(-n_Cr_PAC2.*lambda_Cr_PAC2.*ShearRate.*mu_0_PAC2+n_Cr_PAC2.*lambda_Cr_PAC2.*ShearRate.*mu_inf_PAC2+lambda_Cr_PAC2.*ShearRate.*(lambda_Cr_PAC2.*ShearRate+1).^n_Cr_PAC2.*mu_inf_PAC2+lambda_Cr_PAC2.*ShearRate.*mu_0_PAC2-lambda_Cr_PAC2.*ShearRate.*mu_inf_PAC2+mu_inf_PAC2.*(lambda_Cr_PAC2.*ShearRate+1).^n_Cr_PAC2+mu_0_PAC2-mu_inf_PAC2)./(lambda_Cr_PAC2.*ShearRate.*(lambda_Cr_PAC2.*ShearRate+1).^n_Cr_PAC2.*mu_inf_PAC2+lambda_Cr_PAC2.*ShearRate.*mu_0_PAC2-lambda_Cr_PAC2.*ShearRate.*mu_inf_PAC2+mu_inf_PAC2.*(lambda_Cr_PAC2.*ShearRate+1).^n_Cr_PAC2+mu_0_PAC2-mu_inf_PAC2);
% plot(ShearRate,k_PL.*ShearRate.^(n_PL-1));

% Time scales with local PL coefficients
lambda_FNSD_PAC2 = (PL_FNSD_PAC2.K./(2.*k_PL)).^(1./(PL_FNSD_PAC2.n-n_PL));
% Element-wise n and K with fit function of n and K
% lambda_FNSD_PAC2 = (K_coeff.a1.*exp(-((ShearRate-K_coeff.b1)./K_coeff.c1).^2)./(2.*k_PL)).^(1./(n_coeff.a1.*exp(-((ShearRate-n_coeff.b1)./n_coeff.c1).^2)-n_PL));

FNSD_K_PAC2 = PL_FNSD_PAC2.K;
FNSD_n_PAC2 = PL_FNSD_PAC2.n;

parentpath = cd(cd('..'));
path1 = [parentpath '\Timescales'];
save([path1 '\Timescales_FNSD_PAC2.mat'],...
    'lambda_FNSD_PAC2',...
    'ShearRate',...
    'FNSD_K_PAC2',...
    'FNSD_n_PAC2');

yyaxis right;
g(4) = plot(ShearRate,lambda_FNSD_PAC2,'-','Color',col_PAC2,'LineWidth',2);


%% PAC4
yyaxis left;
% Apparent viscosity estimated, Cox-Merz rule
g(5) = plot(AngularFrequency_PAC4,ComplexViscosity_PAC4,'o','Color',col_PAC4,'MarkerFaceColor',col_PAC4);

% Apparent viscosity, Cross fit
ShearRate = [logspace(-1,3) 1200];
ApparentViscosity=mu_inf_PAC4+(mu_0_PAC4-mu_inf_PAC4)./(1+(lambda_Cr_PAC4.*ShearRate).^n_Cr_PAC4);
g(6) = plot(ShearRate,ApparentViscosity,'-','Color',col_PAC4);

% First normal stress coefficient, Launs rule
FirstNormalStressCoeff_PAC4 = 2.*(StorageModulus_PAC4./(AngularFrequency_PAC4.^2)).*(1+(StorageModulus_PAC4./LossModulus_PAC4).^2).^0.5;
FirstNormalStressDiff_PAC4 = FirstNormalStressCoeff_PAC4.*AngularFrequency_PAC4.^2;
%g(8) = plot(AngularFrequency_PAC4,FirstNormalStressCoeff_PAC4,'d','Color',col_PAC4,'MarkerFaceColor',col_PAC4);

% PL-fit of FNSD & FNSC
PL_FNSD_PAC4 = createFit_PowerLaw(AngularFrequency_PAC4(DataRange_PAC4),FirstNormalStressDiff_PAC4(DataRange_PAC4));
close gcf;
FNSC_fit = PL_FNSD_PAC4.K.*ShearRate.^(PL_FNSD_PAC4.n-2);
% Element-wise n and K with fit function of n and K
% [ FNSC_fit, K_coeff, n_coeff ] = get_FNSC_fit_variable_N_K( AngularFrequency_PAC4(DataRange_PAC4),FirstNormalStressDiff_PAC4(DataRange_PAC4),ShearRate );
g(7) = plot(ShearRate,FNSC_fit,'--','Color',col_PAC4);

% xdata = AngularFrequency_PAC2
% ydata = FirstNormalStressDiff_PAC2
% 
% xdata = AngularFrequency_PAC4
% ydata = FirstNormalStressDiff_PAC4

% Compute local Power Law n & K for current psi & shear rate based on analytical solution
k_PL = exp(log(ShearRate).*(mu_0_PAC4-mu_inf_PAC4).*n_Cr_PAC4.*lambda_Cr_PAC4.*ShearRate./((mu_inf_PAC4.*(lambda_Cr_PAC4.*ShearRate+1).^n_Cr_PAC4+mu_0_PAC4-mu_inf_PAC4).*(lambda_Cr_PAC4.*ShearRate+1))).*(mu_inf_PAC4+(lambda_Cr_PAC4.*ShearRate+1).^(-n_Cr_PAC4)*mu_0_PAC4-(lambda_Cr_PAC4.*ShearRate+1).^(-n_Cr_PAC4).*mu_inf_PAC4);
n_PL=(-n_Cr_PAC4.*lambda_Cr_PAC4.*ShearRate.*mu_0_PAC4+n_Cr_PAC4.*lambda_Cr_PAC4.*ShearRate.*mu_inf_PAC4+lambda_Cr_PAC4.*ShearRate.*(lambda_Cr_PAC4.*ShearRate+1).^n_Cr_PAC4.*mu_inf_PAC4+lambda_Cr_PAC4.*ShearRate.*mu_0_PAC4-lambda_Cr_PAC4.*ShearRate.*mu_inf_PAC4+mu_inf_PAC4.*(lambda_Cr_PAC4.*ShearRate+1).^n_Cr_PAC4+mu_0_PAC4-mu_inf_PAC4)./(lambda_Cr_PAC4.*ShearRate.*(lambda_Cr_PAC4.*ShearRate+1).^n_Cr_PAC4.*mu_inf_PAC4+lambda_Cr_PAC4.*ShearRate.*mu_0_PAC4-lambda_Cr_PAC4.*ShearRate.*mu_inf_PAC4+mu_inf_PAC4.*(lambda_Cr_PAC4.*ShearRate+1).^n_Cr_PAC4+mu_0_PAC4-mu_inf_PAC4);
% plot(ShearRate,k_PL.*ShearRate.^(n_PL-1));

% Time scales with local PL coefficients
lambda_FNSD_PAC4 = (PL_FNSD_PAC4.K./(2.*k_PL)).^(1./(PL_FNSD_PAC4.n-n_PL));
% Element-wise n and K with fit function of n and K
% lambda_FNSD_PAC4 = (K_coeff.a1.*exp(-((ShearRate-K_coeff.b1)./K_coeff.c1).^2)./(2.*k_PL)).^(1./(n_coeff.a1.*exp(-((ShearRate-n_coeff.b1)./n_coeff.c1).^2)-n_PL));


FNSD_K_PAC4 = PL_FNSD_PAC4.K;
FNSD_n_PAC4 = PL_FNSD_PAC4.n;

parentpath = cd(cd('..'));
path1 = [parentpath '\Timescales'];
save([path1 '\Timescales_FNSD_PAC4.mat'],...
    'lambda_FNSD_PAC4',...
    'ShearRate',...
    'FNSD_K_PAC4',...
    'FNSD_n_PAC4');

yyaxis right;
g(8) = plot(ShearRate,lambda_FNSD_PAC4,'-','Color',col_PAC4,'LineWidth',2);

%% Legend
% 
% h_legend = legend( g,...
%     'PAC2 - Complex viscosity $\eta^{*}(\omega = \dot \gamma)$, FS data (Cox-Merz 1958)', ...
%     'PAC2 - Apparent viscosity $\eta (\dot \gamma)$, Cross fit of FC data', ...    
%     ['PAC2 - FNSC $\psi_{1}(\dot \gamma = \omega)$, Power law fit (n = ' num2str(PL_FNSD_PAC2.n,3) ', K = ' num2str(PL_FNSD_PAC2.K,3) ')'],...
%     'PAC2 - Time scale $\lambda_{rs}^{PAC2} = \frac{K_{FNSD}}{2K_{PL}}^\frac{1}{n_{FNSD}-n_{PL}}$ ',...
%     'PAC4 - Complex viscosity $\eta^{*}(\omega = \dot \gamma)$, FS data (Cox-Merz 1958)', ...
%     'PAC4 - Apparent viscosity $\eta (\dot \gamma)$, Cross fit of FC data', ...    
%     ['PAC4 - FNSC $\psi_{1}(\dot \gamma = \omega)$; Power law fit (n = ' num2str(PL_FNSD_PAC4.n,3) ', K = ' num2str(PL_FNSD_PAC4.K,3) ')'],...
%     'PAC4 - Time scale $\lambda_{rs}^{PAC4} = \frac{K_{FNSD}}{2K_{PL}}^\frac{1}{n_{FNSD}-n_{PL}}$ ');

% h_legend = legend( g,...
%     'PAC2 - Complex viscosity $\eta^{*}(\omega = \dot \gamma)$, FS data (Cox-Merz 1958)', ...
%     'PAC2 - Apparent viscosity $\eta (\dot \gamma)$, Cross fit of FC data', ...    
%     'PAC2 - FNSC $\psi_{1}(\dot \gamma = \omega)$, obtained from FS data (Laun 1986)', ...
%     ['PAC2 - FNSD $\psi_{1}(\dot \gamma = \omega)$, Power law fit (n = ' num2str(PL_FNSD_PAC2.n,3) ', K = ' num2str(PL_FNSD_PAC2.K,3) ')'],...
%     'PAC2 - Time scale $\lambda_{rs}^{PAC2} = \frac{K_{FNSD}}{2K_{PL}}^\frac{1}{n_{FNSD}-n_{PL}}$ ',...
%     'PAC4 - Complex viscosity $\eta^{*}(\omega = \dot \gamma)$, FS data (Cox-Merz 1958)', ...
%     'PAC4 - Apparent viscosity $\eta (\dot \gamma)$, Cross fit of FC data', ...    
%     'PAC4 - FNSC $\psi_{1}(\dot \gamma = \omega)$; obtained from FS data (Laun 1986)', ...
%     ['PAC4 - FNSD $\psi_{1}(\dot \gamma = \omega)$; Power law fit (n = ' num2str(PL_FNSD_PAC4.n,3) ', K = ' num2str(PL_FNSD_PAC4.K,3) ')'],...
%     'PAC4 - Time scale $\lambda_{rs}^{PAC4} = \frac{K_{FNSD}}{2K_{PL}}^\frac{1}{n_{FNSD}-n_{PL}}$ ');

% set(h_legend,...
%     'Location','southeast',...
%     'Interpreter','latex');

str = {'\eta^{*}(\omega = \gamma), PAC2',...
'\eta (\gamma), PAC2',...
'\psi_{1}(\gamma = \omega), PAC2',...
'\eta^{*}(\omega = \gamma), PAC4',...
'\eta (\gamma), PAC4',...
'\psi_{1}(\gamma = \omega), PAC4',...
'\lambda_{rs}, PAC2',...
'\lambda_{rs}, PAC4'};
h_legend = columnlegend(2,str);
legend boxoff
pos=get(legend,'position');
set(legend,'position',[0.45*pos(1) 0.25*pos(2) pos(3) pos(4)]);

    
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
set(fig2,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig2.PaperPosition;
fig2.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
path1 = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = 'Merz_NormalStressCoefficient';
print(fig2,[path1 fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig2,[path1 fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig2,[path1 fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig2,[path1 num2str(9)],'-dpng','-r600')
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file
