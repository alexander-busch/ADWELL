function [fitresult, gof] = createFit_Cross(ShearRate, AppViscosity)
%CREATEFIT1(SHEARRATE,APPVISCOSITY)
%  Create a fit.
%
%  Data for 'Cross - Apparent Viscosity' fit:
%      X Input : ShearRate
%      Y Output: AppViscosity
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Apr-2017 11:06:24


%% Fit: 'Cross - Apparent Viscosity'.
[xData, yData] = prepareCurveData( ShearRate, AppViscosity );

% Set up fittype and options.
ft = fittype( 'mu_inf+(mu_0-mu_inf)/(1+(lambda*x)^n)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0.00102 0];
opts.Robust = 'LAR';
opts.StartPoint = [0.05 1 0.00102 0.8];
opts.Upper = [10 10 0.00102 1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Cross - Apparent Viscosity' );
h = plot( fitresult, xData, yData );
legend( h, 'AppViscosity vs. ShearRate', 'Cross - Apparent Viscosity', 'Location', 'NorthEast' );
% Label axes
xlabel ShearRate
ylabel AppViscosity
grid on
set(gca,...
    'XScale','log',...
    'YScale','log');