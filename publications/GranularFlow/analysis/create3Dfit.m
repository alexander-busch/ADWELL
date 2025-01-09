function [fitresult, gof] = create3Dfit(x, y, vof, ft)
%CREATEFIT(X,Y,VOF)
%  Create a fit.
%
%  Data for 'untitled fit 2' fit:
%      X Input : x
%      Y Input : y
%      Z Output: vof
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 11-Jul-2018 17:16:32


%% Fit: 'untitled fit 2'.
[xData, yData, zData] = prepareSurfaceData( x, y, vof );

% Set up fittype and options.
% ft = 'cubicinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'untitled fit 2' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 2', 'vof vs. x, y', 'Location', 'NorthEast' );
% Label axes
xlabel x
ylabel y
zlabel vof
grid on


