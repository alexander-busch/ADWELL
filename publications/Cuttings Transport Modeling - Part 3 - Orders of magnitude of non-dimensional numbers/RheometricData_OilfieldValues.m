%% Shear stress and apparent viscosity based on Oilfield Bingham model and Fann viscometer

% Representative oilfield PV & YP acc. to Busch et al. (2016)
YP = [0; 0; 10/0.4788026; 10/0.4788026]; % [Pa] --> [lbf/100ft²]
PV = [15; 30; 15; 30]; % [mPa.s]

% Shear rate range
SR = logspace(-2,4);

% Plot
fig = figure('Name','Re = f(gpm)','color','w');
hold on; grid on; box on;

for ii=1:4
    
    tau = 0.4788026 * 1.066.*(YP(ii)+PV(ii)./(1022-511).*SR);
    eta = tau./SR;
    plot(SR,eta);
end

xlabel 'Shear rate [1/s]';
ylabel 'Apparent viscosity [Pa.s]';
set(gca,...
    'XScale','log',...
    'YScale','log');