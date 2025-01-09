%% Import data

% Shear rate
ShearRate = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',2,'A1:A50');

% Shear stress
Tau1 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',2,'B1:B50');
Tau2 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',3,'B1:B50');
Tau3 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',4,'B1:B50');

% App. viscosity
Vis1 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',2,'C1:C50');
Vis2 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',3,'C1:C50');
Vis3 = xlsread('M:\data\AdWell\1 - Rheology\Rheology model\Glydril.xlsx',4,'C1:C50');


%% Analyze data

% Shear stress
Tau=ones(3,length(Tau1));
Tau(1,:)=Tau1;
Tau(2,:)=Tau2;
Tau(3,:)=Tau3;
Tau_M  = mean(Tau,1); % Apparent viscosity mean
Tau_SD = std(Tau,1); % Apparent viscosity standard deviation

% App. viscosity
Vis=ones(3,length(Vis1));
Vis(1,:)=Vis1;
Vis(2,:)=Vis2;
Vis(3,:)=Vis3;
Vis_M  = mean(Vis,1); % Apparent viscosity mean
Vis_SD = std(Vis,1); % Apparent viscosity standard deviation


%% Plot data


figure('name','Shear stress vs. Shear rate');
hold on;
plot(ShearRate, Tau_M,'g','LineWidth',1);
plot(ShearRate, Tau_M + 3*Tau_SD,'g--');
plot(ShearRate, Tau_M - 3*Tau_SD,'g--');
set(gca,'XScale','log','YScale','log');

figure('name','App. viscosity vs. Shear rate');
hold on;
plot(ShearRate, Vis_M,'g','LineWidth',1);
plot(ShearRate, Vis_M + 3*Vis_SD,'g--');
plot(ShearRate, Vis_M - 3*Vis_SD,'g--');
set(gca,'XScale','log','YScale','log');

%% Plot fit

ny_0 = 169.6
ny_inf = 3.827e-2
K = 0.01329
k = 1
n = 0.6203

Fit=(ny_0*K+ny_inf.*(k.*ShearRate).^n)./(K+(k.*ShearRate).^n);
plot(ShearRate, Fit);
