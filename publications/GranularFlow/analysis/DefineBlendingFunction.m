close all;
clear all;
clc;
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';

fig= CreateFigure( '', 'Solid volume fraction \alpha_s [-]', 'Frictional pressure p_{s,fric} [Pa]', {'lin','log'}, 'DINA4' );
TightFigure(gca);

alpha_s = linspace(0.55,0.63);
loc = 0.6;
dis = 0.001;

f1=0.1*alpha_s.*(alpha_s-min(alpha_s)).^2./((max(alpha_s)-alpha_s).^5);
f2=alpha_s.*10.^25.*(alpha_s-min(alpha_s)).^10;



plot(alpha_s,f1);
plot(alpha_s,f2);

yyaxis right;

blend = (1+ tanh((alpha_s-loc)/dis))/2;
%blend = tanh((alpha_s-loc)/dis)/2;
plot(alpha_s,blend);


yyaxis left;
result = (1 - blend).*f1 + (blend.*f2);
plot(alpha_s,result);