function [fitresult, gof] = createFit_CarreauYasuda(ShearRate, AppViscosity)
%CREATEFIT1(SHEARRATE,APPVISCOSITY)
%  Create a fit.
%
%  Data for 'Carreau-Yasuda' fit:
%      X Input : ShearRate
%      Y Output: AppViscosity
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Apr-2017 11:13:17


%% Fit: 'Carreau-Yasuda'.
[xData, yData] = prepareCurveData( ShearRate, AppViscosity );

% Set up fittype and options.
ft = fittype( 'ny_inf+(ny_0-ny_inf)*(1+(lambda*x)^a)^((n-1)/a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 0 0 0.00102];
opts.StartPoint = [0.791815130930056 0.823457828327293 0.9502 0.694828622975817 0.00102];
opts.Upper = [Inf 100 1 10 0.00102];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Carreau-Yasuda' );
h = plot( fitresult, xData, yData );
legend( h, 'AppViscosity vs. ShearRate', 'Carreau-Yasuda', 'Location', 'NorthEast' );
% Label axes
xlabel ShearRate
ylabel AppViscosity
grid on
set(gca,...
    'XScale','log',...
    'YScale','log');


