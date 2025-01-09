function [ realtmax ] = GetRealFlowTime( tmax )
%tmax = 30;
d_o = 0.216; % Outer diameter
d_j = 0.168; % Tool joint diameter
rpm_pipe = 130;
rpm_whirl = d_j/d_o*rpm_pipe;
dt = 1e-3;

% Duration of one orbital cycle
T = round(1/(rpm_whirl/60),3);

% Replace mesh after n orbital cycles
n = 4;

% Replace mesh after time_replacemesh seconds
time_replacemesh = round(n*T,3);

% Set timestep
timestep = round(time_replacemesh/dt,0);

% Number of required inner loops to achieve total flow time
jMax =  round(tmax/time_replacemesh,0);

% Real flowtime
realtmax = jMax*time_replacemesh;
realtmax = (jMax-1)*time_replacemesh;
end

