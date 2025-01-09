%% Clean-up
clear all;
close all;
clc;

addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';

% Parameters
global rho_f rho_s lambda_Cross n_Cross mu_0 mu_inf FNSC_n FNSC_K;
global lambda_el_upper lambda_el_lower lambda_th_upper lambda_th_lower;
global fig_MFS_3ITT fig_MFS_Maxwell fig_PS_3ITT fig_PS_Maxwell fig_MFS_Pipkin_th fig_MFS_Pipkin_el fig_PS_Pipkin_th fig_PS_Pipkin_el colorlist;

fluid = 'PAC4';

% Outer & inner diameter [inch] combinations acc. to K&M
d_o = [17.5; 12.25; 9.875; 8.5; 6.125];
d_i = [6.625; 6.625; 5; 4.5; 3.5];

% Hydraulic diameter [meter]
d_h = (d_o - d_i) * 25.4/1000;

% Length scale of Main Flow Scale
L = 1; % [1; 10; 100; 1000];

% Volumetric flow rate as f(d_h) [lpm] --> [m³/s] acc. to K&M
n = 10;
VolFlowRate = logspace(1,3,n)/1000/60;% Volumetric flow rate range [lpm] --> [m³/s]
VolFlowRateMax = [5700; 4200; 3500; 2300; 700] / 1000/60;

% Densities
rho_f = 1000;
rho_s = 2650;

% Rotational speed of drill pipe
rpm = [50; 100; 150; 200];

% Particle diameter [mm] --> [m]
d_p = [0.01; 0.1; 1; 10] /1000;

% Settling distance [m]
H = max(d_h);

% Stokes settling velocity guess (m/s)
v_set = [4e-7; 4e-5; 4e-3; 7e-1];

% Get fluid data
if strcmp(fluid, 'PAC2')
    % Load rheological time scales
    load FC_Cross_PAC2.mat;
    load Timescales_FC_PAC2.mat;
    load Timescales_FS.mat;
    load Timescales_3ITT.mat;
    load Timescales_3ITT_Large.mat;
    load Timescales_FNSD_PAC2.mat;
    % Assign parameters
    lambda_el_lower = min(lambda_Cr_PAC2, min(MaxwellRelaxationTime_PAC2));
    lambda_el_upper = max(lambda_Ca_PAC2, max(MaxwellRelaxationTime_PAC2));
    lambda_th_lower = min(lambda_3ITT_100_50_PAC2, lambda_3ITT_300_600_PAC2);
    lambda_th_upper = max(lambda_3ITT_100_50_PAC2, lambda_3ITT_300_600_PAC2);
    lambda_Cross = lambda_Cr_PAC2;
    n_Cross = n_Cr_PAC2;
    mu_0 = mu_0_PAC2;
    mu_inf = mu_inf_PAC2;
    FNSC_n = FNSD_n_PAC2;
    FNSC_K = FNSD_K_PAC2;
else
    % Load rheological time scales
    load FC_Cross_PAC4.mat;
    load Timescales_FC_PAC4.mat;
    load Timescales_FS.mat;
    load Timescales_3ITT.mat;
    load Timescales_3ITT_Large.mat;
    load Timescales_FNSD_PAC4.mat;
    % Assign parameters
    lambda_el_lower = min(lambda_Cr_PAC4, min(MaxwellRelaxationTime_PAC4));
    lambda_el_upper = max(lambda_Ca_PAC4, max(MaxwellRelaxationTime_PAC4));
    lambda_th = [lambda_3ITT_100_50_PAC4, lambda_3ITT_100_50_PAC4_LowShearRate, lambda_3ITT_300_600_PAC4, lambda_3ITT_300_600_PAC4_HighShearRate, lambda_3ITT_300_600_PAC4_LowShearRate];
    lambda_th_lower = min(lambda_th);
    lambda_th_upper = max(lambda_th);
    lambda_Cross = lambda_Cr_PAC4;
    n_Cross = n_Cr_PAC4;
    mu_0 = mu_0_PAC4;
    mu_inf = mu_inf_PAC4;
    FNSC_n = FNSD_n_PAC4;
    FNSC_K = FNSD_K_PAC4;
end


%% Create figures

% Colors
% http://se.mathworks.com/help/matlab/ref/colorspec.html
colorlist = {[0.3 0.3 0.3] ...
    [0.6 0.6 0.6] ...
    [0.9 0.9 0.9] ... % Blue [0 0 1]
    'k'}; 

colorlist = {[0.0 0.0 0.0] ...
    [0.4 0.4 0.4] ...
    [0.8 0.8 0.8] ... % Blue [0 0 1]
    'k'}; 


% Main Flow Scale - 3ITT
fig_MFS_3ITT  = figure('color','w','Name','Degree of Thixotropy, Annulus scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e0 1e4],...
    'ylim', [1e-3 1e5],...
    'box','on');

xlabel('Re_G [-]');
ylabel('De_t_h, Wi_t_h [-]');


% Main Flow Scale - Maxwell
fig_MFS_Maxwell = figure('color','w','Name','Degree of Viscoelasticity, Annulus scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e0 1e4],...
    'ylim', [1e-4 1e2],...
    'box','on');
    % ,...

xlabel('Re_G [-]');
ylabel('De_e_l, Wi_e_l [-]');

% Particle Scale - 3ITT
fig_PS_3ITT = figure('color','w','Name','Degree of Thixotropy, Particle scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-7 1e3],...
    'ylim', [1e-4 1e4],...
    'box','on');
    % ,...

xlabel('Re_p [-]');
ylabel('De_t_h, Wi_t_h [-]');

% Particle Scale - Maxwell
fig_PS_Maxwell = figure('color','w','Name','Degree of Viscoelasticity, Particle scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-7 1e3],...
    'ylim', [1e-4 1e4],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('Re_p [-]');
ylabel('De_e_l, Wi_e_l [-]');


% Main Flow Scale - Pipkin
fig_MFS_Pipkin_th = figure('color','w','Name','Pipkin space, Main flow scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-3 1e3],...
    'ylim', [1e-1 1e5],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('De_t_h [-]');
ylabel('Wi_t_h [-]');

fig_MFS_Pipkin_el = figure('color','w','Name','Pipkin space, Main flow scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-5 1e0],...
    'ylim', [1e-3 1e2],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('De_e_l [-]');
ylabel('Wi_e_l [-]');


% Particle Scale - Pipkin
fig_PS_Pipkin_th = figure('color','w','Name','Pipkin space, Particle scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [2e-1 6e3],...
    'ylim', [4e-1 2e4],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('De_t_h [-]');
ylabel('Wi_t_h [-]');

fig_PS_Pipkin_el = figure('color','w','Name','Pipkin space, Particle scale','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;

set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [6e-3 2e1],...
    'ylim', [1e-2 3e1],...
    'box','on');
    % ,...
    % 'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    % 'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    % 'FontSize',24);

xlabel('De_e_l [-]');
ylabel('Wi_e_l [-]');

%% Evaluate time scales & De

for ii = 1:length(d_h)
    
    % Assign max. vol. flow rate for current d_h to vol. flow rate range
    VolFlowRate(n) = VolFlowRateMax(ii);
    
    % Evaluate timescales & plot De vs. Re
    CompMainFlowScale(d_o(ii), d_i(ii), d_h(ii), L, VolFlowRate, fluid, ii, 1 );
   
    % Debugging
    % d_o = d_o(ii)
    % d_i = d_i(ii)
    % d_h = d_h(ii)
end


% for i = 1:length(d_p)
   
    % Evaluate timescales & plot De vs. Re
    CompParticleScale( d_p, fluid, v_set );
    
    % Debugging
    % d_p = d_p(i)
    % v_set = v_set(i)
%end



%% Format figures
fig_path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2018-05 - Judging the GNF assumption for cuttings transport modeling by applying time scale comparions\Figures\';

% Main Flow Scale - 3ITT
figure(fig_MFS_3ITT);
% legend('Mean Fligth time',...
%     'Large Eddy Turnover time',...
%     'Shear Rate',...
%     'Rotation',...
%     'Kolmogorov',...
%     'location','northwest');
legend('De_t_h',...
    'Wi_t_h',...
    'location','northwest');

hline(1e1,'k-');
hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_MFS_3ITT'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Main Flow Scale - Maxwell
figure(fig_MFS_Maxwell);
% legend('Mean flight time',...
%     'Large Eddy Turnover time',...
%     'Shear Rate',...
%     'Rotation',...
%     'Kolmogorov',...
%     'location','northwest');
legend('De_e_l',...
    'Wi_e_l',...
    'location','northwest');
    hline(1e1,'k-');
    hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig.PaperUnits = 'centimeters';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_MFS_Maxwell'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Particle Scale - 3ITT
figure(fig_PS_3ITT);
legend('De_t_h',...
    'We_t_h',...
    '\tau_t_h/T_p_,_s_e_t',...
    'location','northwest');
    hline(1e1,'k-');
    hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_PS_3ITT'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Particle Scale - Maxwell
figure(fig_PS_Maxwell);
legend('De_e_l',...
    'We_e_l',...
    '\tau_e_l/T_p_,_s_e_t',...
    'location','northwest');
    hline(1e1,'k-');
    hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_PS_Maxwell'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Main Flow Scale - Pipkin - th
figure(fig_MFS_Pipkin_th);
% legend('Settling time',...
%     'Stokes relaxation time',...
%     'Shear Rate',...
%     'location','northwest');
%     hline(1e1,'k-');
%     hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_MFS_Pipkin_th'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Main Flow Scale - Pipkin - el
figure(fig_MFS_Pipkin_el);
% legend('Settling time',...
%     'Stokes relaxation time',...
%     'Shear Rate',...
%     'location','northwest');
%     hline(1e1,'k-');
%     hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = [fluid '_MFS_Pipkin_el'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Particle Scale - Pipkin - th
figure(fig_PS_Pipkin_th);
% legend('Settling time',...
%     'Stokes relaxation time',...
%     'Shear Rate',...
%     'location','northwest');
%     hline(1e1,'k-');
%     hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_name = [fluid '_PS_Pipkin_th'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Main Flow Scale - Pipkin - el
figure(fig_PS_Pipkin_el);
% legend('Settling time',...
%     'Stokes relaxation time',...
%     'Shear Rate',...
%     'location','northwest');
%     hline(1e1,'k-');
%     hline(1e-1,'k-');
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
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = [fluid '_PS_Pipkin_el'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file