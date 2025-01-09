function [ fig ] = CreateFigure( fig_title, label_x, label_y )
%	CreateFigure Create and format figure 
%   Detailed explanation goes here

% fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 21.0 14.8]); % DIN A5 ;
fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 29.7 21.0]); % DIN A4 ;
% fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 42.0 29.7]); % DIN A3 ;
% fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 59.4 42.0]); % DIN A2 ;
% fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 84.1 59.4]); % DIN A1 ;
% fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',[1 1 118.9 84.1]); % DIN A0 ;

hold on;
grid on;
box on;

% title(fig_title);
xlabel(label_x);
ylabel(label_y);

% set(gca,...
%     'XScale','log',...
%     'YScale','log');



end

