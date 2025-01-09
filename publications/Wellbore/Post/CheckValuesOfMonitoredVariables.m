close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));

filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\remote';
%filepath_simdata = 'C:\Users\alexabus\Desktop\hpcresults';
filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive\181027_1\5';

% fullpath = [filepath_simdata '\solution\monitors\fluid-x-mass-flow-rate.txt'];
% [time, m_dot_f] = importMonitoredQuantities(fullpath);
% 
% fullpath = [filepath_simdata '\solution\monitors\solid-x-mass-flow-rate.txt'];
% [time, m_dot_s] = importMonitoredQuantities(fullpath);

fullpath = [filepath_simdata '\solution\monitors\fluid-superficial-x-velocity.txt'];
[time, U_f] = importMonitoredQuantities(fullpath);

fullpath = [filepath_simdata '\solution\monitors\solid-superficial-x-velocity.txt'];
[time, U_s] = importMonitoredQuantities(fullpath);

dpdx = (100:100:500)';

dt = 1e-3;
tmax = GetRealFlowTime( 30 );
idx_low = 20*(1/dt);
idx_up = tmax*(1/dt);


%% Plot vs time
CreateFigure('U_f, U_s, CTR vs. time','Time [s]','U_f, U_s  [m/s]', {'lin' 'lin'}, 'DINA5');

yyaxis left
% plot(time,m_dot_f)
% plot(time,m_dot_s)
% set(gca,...
%     'xlim', [0 round(time(end),0)],...
%     'ylim', [0 round(m_dot_f(end),0)]);
plot(time,U_f)
plot(time,U_s)
set(gca,...
    'xlim', [0 round(time(end),0)]);
% ,...
%     'ylim', [0 round(U_f(end),2)]

yyaxis right
time2=dt*[0; repelem((1:1:length(dpdx)-1)',2,1); length(dpdx)]*idx_up;
dpdx2 = repelem(dpdx,2,1);
plot(time2,dpdx2);
set(gca,...
    'xlim', [0 round(time(end),0)],...
    'ylim', [0 dpdx2(end)]);
ylabel('-\Deltap/\Deltax [Pa/m]');

TightFigure(2);


%% Plot vs. dpdx

mean_U_f = zeros(length(dpdx),1);
mean_U_s = mean_U_f; CTR = mean_U_f;

for ii = 1:length(dpdx)
    mean_U_f(ii) = mean(U_f(ii*idx_low:ii*idx_up));
    mean_U_s(ii) = mean(U_s(ii*idx_low:ii*idx_up));
    CTR(ii) = mean_U_s(ii)/(mean_U_f(ii)+mean_U_s(ii));
    CTR(ii) = mean_U_s(ii)/mean_U_f(ii);
end

CreateFigure('U_f, U_s, CTR vs. dpdx','Pressure gradient -\Deltap/\Deltax [Pa/m]','Time averaged U_f, U_s  [m/s]', {'lin' 'lin'}, 'DINA5');

yyaxis left
plot(dpdx,mean_U_f);
plot(dpdx,mean_U_s);
set(gca,...
    'xlim', [50 550],...
    'ylim', [0 2]);

yyaxis right
plot(dpdx,CTR)
set(gca,...
    'xlim', [50 550],...
    'ylim', [0 1]);
ylabel('CTR [-]');

TightFigure(2);

%%

% filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\hpc';
% 
% fullpath = [filepath_simdata '\solution\monitors\fluid-x-mass-flow-rate.txt'];
% [time, m_dot_f] = importMonitoredQuantities(fullpath);
% 
% fullpath = [filepath_simdata '\solution\monitors\solid-x-mass-flow-rate.txt'];
% [time, m_dot_s] = importMonitoredQuantities(fullpath);
% 
% fullpath = [filepath_simdata '\solution\monitors\fluid-superficial-x-velocity.txt'];
% [time, U_f] = importMonitoredQuantities(fullpath);
% 
% fullpath = [filepath_simdata '\solution\monitors\solid-superficial-x-velocity.txt'];
% [time, U_s] = importMonitoredQuantities(fullpath);
% 
% yyaxis left
% plot(time,m_dot_f,':','LineWidth',2)
% plot(time,m_dot_s,':','LineWidth',2)
% set(gca,...
%     'xlim', [0 round(time(end),0)],...
%     'ylim', [0 round(m_dot_f(end),0)]);
% 
% yyaxis right
% hold on
% plot(time,U_f,':','LineWidth',2)
% plot(time,U_s,':','LineWidth',2)
% set(gca,...
%     'xlim', [0 round(time(end),0)],...
%     'ylim', [0 round(U_f(end),0)]);

