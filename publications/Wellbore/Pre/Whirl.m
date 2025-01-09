% Whirl

clear all;
close all;
clc;

% Parameter
d_o = 0.04;
d_i = 0.025;
rpm = [0 100 200 300];
Q_f = [0; 0.00026; 0.00041; 0.00094];
[X, Y] = meshgrid(rpm,Q_f);

% Experimental data from Milad Khatibi ("Position of Drill String.xlsx")
% rows = f(rpm)
% columns = f(Q_f)
eccentricity = [-0.0071 -0.0065 -0.0058 -0.0045;...
    NaN -0.0064 -0.0057 -0.0042;
    NaN -0.0064 -0.0050 NaN;...
    NaN -0.0056 NaN NaN];

amplitude = [0.0000 0.0046 0.0038 0.0016;...
    NaN 0.0041 0.0043 0.0022;
    NaN 0.0038 0.0027 NaN;...
    NaN 0.0035 NaN NaN];

% Assumption: For rpm = 0, the datapoints for Q_f show the same trend as
% for Q_f = 0
eccentricity(:,1) = eccentricity(1,1)+abs(eccentricity(1,2)-eccentricity(:,2));
amplitude(:,1) = 0;


% Plot
figure; 
hold on;

scatter3(reshape(X,[numel(X),1]), reshape(Y,[numel(X),1]), reshape(eccentricity,[numel(X),1]));
scatter3(reshape(X,[numel(X),1]), reshape(Y,[numel(X),1]), reshape(amplitude,[numel(X),1]));
scatter3(reshape(X,[numel(X),1]), reshape(Y,[numel(X),1]), reshape(eccentricity + amplitude,[numel(X),1]));

[coeff_eccentricity, gof] = createPolynomialFit(rpm, Q_f, eccentricity);
fit_eccentricity = coeff_eccentricity.p00 + coeff_eccentricity.p10.*X + coeff_eccentricity.p01.*Y + coeff_eccentricity.p20.*X.^2 + coeff_eccentricity.p11.*X.*Y + coeff_eccentricity.p02.*Y.^2;
surf(X, Y, fit_eccentricity);

[coeff_amplitude, gof] = createPolynomialFit(rpm, Q_f, amplitude);
fit_amplitude = coeff_amplitude.p00 + coeff_amplitude.p10.*X + coeff_amplitude.p01.*Y + coeff_amplitude.p20.*X.^2 + coeff_amplitude.p11.*X.*Y + coeff_amplitude.p02.*Y.^2;
surf(X, Y, fit_amplitude);

surf(X, Y, fit_eccentricity+fit_amplitude);

view(40,10)


% Compute individual point
X = 300; % rpm
Y = 0.000251; % Q_f
eccentricity = coeff_eccentricity.p00 + coeff_eccentricity.p10.*X + coeff_eccentricity.p01.*Y + coeff_eccentricity.p20.*X.^2 + coeff_eccentricity.p11.*X.*Y + coeff_eccentricity.p02.*Y.^2
amplitude = coeff_amplitude.p00 + coeff_amplitude.p10.*X + coeff_amplitude.p01.*Y + coeff_amplitude.p20.*X.^2 + coeff_amplitude.p11.*X.*Y + coeff_amplitude.p02.*Y.^2
scatter3(X,Y,eccentricity);
scatter3(X,Y,amplitude);
