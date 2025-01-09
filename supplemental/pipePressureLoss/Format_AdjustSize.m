function [  ] = Format_AdjustSize (  )
%Format_AdjustSize Summary of this function goes here
%   Detailed explanation goes here

% % Adjust size to screen size
% screen_size = get(0, 'ScreenSize');
% origSize = get(gcf, 'Position'); % grab original on screen size
% set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size

% Adjust size to screen size
screen_size = get(0, 'ScreenSize');
origSize = get(gcf, 'Position'); % grab original on screen size
offset_hor = 10;
offset_ver = 80;
set(gcf, 'Position', [offset_hor 0 screen_size(3)-2*offset_hor screen_size(4)-offset_ver ] ); %set to scren size


% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Re-Adjust size
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
% set(fig,'Units','centimeters','Position',[0 0 32 18])
% set(fig,'Position', origSize) %set back to original dimensions


end

