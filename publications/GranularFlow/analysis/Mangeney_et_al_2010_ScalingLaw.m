% Mangeney et al. (2010) scaling law for final deposit height and run-out
% distance

% Scaling law coefficients
k = 0.5; % 0.4
beta = 0.76;
%n_x = 1.0; % For AR < 1.6 acc. to Delannay et al. (2017)
n_y = 2/3;

% Dimensional run-out distance
Cases.x_f_Exp(aa) = Cases.y_0(aa)*2*k/tan(Cases.theta_s(aa)*pi/180); % Equation (5)

% Dimensional deposit height
Cases.y_f_Exp(aa) = Cases.y_0(aa)*beta*(tan(Cases.theta_s(aa)*pi/180)/Cases.AR(aa)/k)^n_y; % Equation (6)

% Plot on dimensional shape figure
figure(fig_Shape);
scatter([min(xlim) Cases.x_f_Exp(aa)],[ Cases.y_f_Exp(aa) 0],'r+');

% Normalized run-out distance
x_n = (Cases.x_f_Exp(aa)-Cases.x_0(aa))/Cases.x_0(aa); % Coordinate system with x0 = 0
x_n = Cases.x_f_Exp(aa)/Cases.x_0(aa); % Coordinate system with x0 = x_0

% Normalized deposit height
y_n = Cases.y_f_Exp(aa)/Cases.x_0(aa);

% Plot on non-dimensional f(AR) figure
figure(fig_AR);
scatter(Cases.AR(aa),x_n,'b');
scatter(Cases.AR(aa),y_n,'g');

% Plot on non-dimensional shape figure
figure(fig_Shape_y_0);
scatter([-1 x_n*x_0/y_0-1],[y_n*x_0/y_0 0],'r');