close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
addpath([pwd '\utilities']);
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Wellbore\Figures\';

CreateFigure( 'Mesh dependency',...
    '', 'Mass flow rate [kg/s]',...
    {'lin' 'lin'},'DINA5');


xtitles = {'Coarse (High Re)' 'Intermediate' 'Fine'	'Superfine (Low Re)'}


yyaxis left
% CFD
cfd=[11.012 11.542 11.414 11.121]
%plotcolor = 'blue'
%plot([1 2 3 4],cfd,'-o','LineWidth',2,'MarkerEdgeColor',plotcolor,'MarkerFaceColor',plotcolor,'Color',plotcolor);
plot([1 2 3 4],cfd,'-o','LineWidth',1);

% Blasius
dpdx = 30
d_o = 0.216
d_i = 0.127
d_h = d_o-d_i
rho=998
eta=0.00102
u=(16*dpdx^4*d_h^5/0.316^4/rho^3/eta)^(1/7)% Blasius (1912)
mdot=u*0.25*pi*(d_o^2-d_i^2)*rho
plot([1 2 3 4],[mdot mdot mdot mdot],'--o','LineWidth',1);

% Haaland
% Re = u*rho*d_h/eta
% f = 1./(4.* (-1.8.*log10((0.00003./3.7/d_h).^1.11+(6.9./Re))).^2);
% f = 1./(1.* (-1.8.*log10((0.00003./3.7/d_h).^1.11+(6.9./Re))).^2);
% dpdx = f.*rho.*u.^2./d_h./2;

set(gca,...
    'FontSize',12,...
    'yLim',[11 12]);


yyaxis right
plot([1 2 3 4],cfd/mdot,'--o','LineWidth',1);
plot([1 2 3 4],cfd/cfd(4),'-o','LineWidth',1); % (cfd(1)+cfd(4))

set(gca,...
    'FontSize',12,...
    'yLim',[0.9 1.1]);



%xlabel(label_x{2},'interpreter','tex','FontName','Arial');
%ylabel(label_y,'interpreter','tex','FontName','Arial');




legend('CFD results','Blasius (1912)','CFD/Blasius','CFD/highRe)');
xticks([1 2 3 4])
xticklabels(xtitles);
TightFigure(gca); % Handles to axis, if only one just provide one
legend('location','south')

% Print to files
fig_name = 'mesh';
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
set( gcf,'PaperType','A5','PaperOrientation','Landscape')
print(gcf,[fig_path fig_name],'-dpdf');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 
