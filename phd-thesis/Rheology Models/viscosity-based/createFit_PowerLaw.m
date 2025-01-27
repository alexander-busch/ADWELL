function [fitresult, gof] = createFit_PowerLaw(ShearRate, AppViscosity)
%CREATEFIT(SHEARRATE,APPVISCOSITY)
%  Create a fit.
%
%  Data for 'Power-Law - Apparent Viscosity' fit:
%      X Input : ShearRate
%      Y Output: AppViscosity
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Apr-2017 10:46:58


%% Fit: 'Power-Law - Apparent Viscosity'.
[xData, yData] = prepareCurveData( ShearRate, AppViscosity );

% Set up fittype and options.
ft = fittype( 'K*x^(n-1)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.546881519204984 0.957506835434298];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Power-Law - Apparent Viscosity' );
h = plot( fitresult, xData, yData );
legend( h, 'AppViscosity vs. ShearRate', 'Power-Law - Apparent Viscosity', 'Location', 'NorthEast' );
% Label axes
xlabel ShearRate
ylabel AppViscosity
grid on
set(gca,...
    'XScale','log',...
    'YScale','log');

