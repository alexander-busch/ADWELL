function [fitresult, gof] = createFit1(xdata, ydata)
%CREATEFIT1(XDATA,YDATA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xdata
%      Y Output: ydata
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Feb-2018 12:36:53


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xdata, ydata );

% Set up fittype and options.
ft = fittype( 'tau_0+K*x^n', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 1 0];
opts.StartPoint = [0.957166948242946 0.8002804688888 0.421761282626275];
opts.Upper = [Inf 3 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'ydata vs. xdata', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel xdata
ylabel ydata
grid on
set(gca,...
    'XScale','log',...
    'YScale','log');
