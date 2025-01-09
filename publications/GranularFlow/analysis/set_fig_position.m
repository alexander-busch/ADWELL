function [  ] = set_fig_position( fighandle, index )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Get PC's active screen size
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

% Matlab does not consider the height of the figure's toolbar...
% or the width of the border... it only cares about the content
toolbar_height = 77;
window_border  = 5;

% Set figure units to pixels
set(fighandle, 'Units', 'pixels');

% get position of figure
figposition = get(fighandle, 'Position');

% Define position
if index==1
    xpos = window_border+0;
    ypos = screen_size(4)-toolbar_height-figposition(4);
elseif index==2
    xpos = window_border+0;
    ypos = 0;
elseif index==3
    xpos = window_border+figposition(3)+window_border;
    ypos = 0;
else
    xpos = 2*(window_border+figposition(3)+window_border);
    ypos = 0;
end

% Reset position of figure [x y w h]
set(fighandle, 'Position',...
    [xpos,...
    ypos,...
    figposition(3),...
    figposition(4)]);
% set(fighandle, 'Position', figposition)

% Reset figure units to centimeters
set(fighandle, 'Units', 'centimeters');


end

