% Orders of magnitude estimate of non-dimensional numbers
% characterizing cuttings transport

close all;
clear all;
clc;


%% Cuttings bed present


fig= CreateFigure( '', 'Slope [°]', 'Force per unit volume f_i [N/m³]', {'lin','lin'}, 'DINA4' );
TightFigure(gca);

phi_rep = 33*pi/180;
phi = sort([linspace(0,90)*pi/180 phi_rep]);

d_s = 0.0001;
rho_s = 2650;
g_mag = 9.81;

V_s = pi*d_s^3/6;


% Gravitational force
plot(phi*180/pi,ones(1,length(phi))*rho_s*g_mag,'k--');

% Frictional force
f_fric = tan(phi_rep)*rho_s*g_mag*cos(phi);
COI = get(gca,'ColorOrderIndex');
plothandle = plot(phi*180/pi,f_fric,'--');
set(get(get(plothandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'ColorOrderIndex',COI);
plot(phi(phi>=phi_rep)*180/pi,f_fric(phi>=phi_rep),'-','linewidth',2);

% Tangential force
f_tan = rho_s*g_mag*sin(phi);
COI = get(gca,'ColorOrderIndex');
plothandle = plot(phi*180/pi,f_tan,'--');
set(get(get(plothandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'ColorOrderIndex',COI);
plot(phi(phi<=phi_rep)*180/pi,f_tan(phi<=phi_rep),'-','linewidth',2);

%x = phi(phi<=phi_rep)*180/pi;
%fill([x fliplr(x)],[f_tan(phi<=phi_rep) fliplr(f_fric(phi<=phi_rep))],[0.8 0.8 0.8])

set(gca,'Xdir', 'reverse');
legend('\rho g',...
    'f_f_r_i_c = tan(\phi_r_e_p) \rho g cos(\phi)',...
    'f_t_a_n = \rho g sin(\phi)',...
    'Location','northeast');





% textstring = ['For respective min. and max. bulk flow rates, recommended drill pipe rotation rates and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2) ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,max(ylim)+200,textstring);

