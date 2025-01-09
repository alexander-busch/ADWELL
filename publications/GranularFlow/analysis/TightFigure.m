function [ ] = TightFigure(ax)
% Expand Axes to Fill Figure (--> Minimum white space)

if length(ax) == 1  
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    factor_hor = 0.015;
    factor_ver = 0.015;
    left = outerpos(1) + ti(1) + factor_hor;
    bottom = outerpos(2) + ti(2) + factor_ver;
    ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
    ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
    ax.Position = [left bottom ax_width ax_height];
else
    outerpos1 = ax(1).OuterPosition; % [left bottom width height]
    outerpos2 = ax(2).OuterPosition; % [left bottom width height]
    ti1 = ax(1).TightInset; % [left bottom right top]
    ti2 = ax(2).TightInset; % [left bottom right top]
    factor_hor = 0.01;
    factor_ver = 0.01;
    left = outerpos1(1) + ti1(1) + factor_hor;
    bottom = outerpos1(2) + ti1(2) + factor_ver;
    ax1_width = outerpos1(3) - ti1(1) - ti2(3) - 2*factor_hor;
    ax1_height = outerpos1(4) - ti1(2) - ti2(4) - 2*factor_ver;
    ax(1).Position = [left bottom ax1_width ax1_height];
    ax(2).Position = [left bottom ax1_width ax1_height];
end
