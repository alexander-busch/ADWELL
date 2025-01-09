function [ fig ] = CreateFigure( fig_title, label_x, label_y )
%	CreateFigure Create and format figure 
%   Detailed explanation goes here

fig = figure('Name',fig_title,'color','w'); % ('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
box on;

% title(fig_title);
xlabel(label_x);
ylabel(label_y);

set(gca,...
    'XScale','log',...
    'YScale','log');

end

