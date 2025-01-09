% close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));



CreateFigure('U_f, U_s vs. time','Time [s]','U_f, U_s [m/s]', {'lin' 'lin'}, 'DINA5');

set(gca,...
    'FontSize',12,...
    'xLim',[0 48],...
	'yLim',[0 3.5]);

pointsintime = [11 16; 27 32; 43 48];

for ii = 1:size(pointsintime,1)
    plot([pointsintime(ii,1) pointsintime(ii,1)],ylim,'k-','HandleVisibility','off')
    plot([pointsintime(ii,2) pointsintime(ii,2)],ylim,'k-','HandleVisibility','off')
end


% AdWell coarse mesh
filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive\181019_3';

fullpath = [filepath_simdata '\solution\monitors\fluid-superficial-x-velocity.txt'];
[time, U_f] = importMonitoredQuantities(fullpath);

fullpath = [filepath_simdata '\solution\monitors\solid-superficial-x-velocity.txt'];
[time, U_s] = importMonitoredQuantities(fullpath);

plot(time,U_f,'--','Color',[0 112/256 192/256],'LineWidth',1);
plot(time,U_s,'--','Color',[204/256 102/256 0],'LineWidth',1);



for ii = 1:size(pointsintime,1)
    mean_U_f_coarse(ii) = mean(U_f(find(abs(time-pointsintime(ii,1)) < 0.0001):find(abs(time-pointsintime(ii,2)) < 0.0001)));
    mean_U_s_coarse(ii) = mean(U_s(find(abs(time-pointsintime(ii,1)) < 0.0001):find(abs(time-pointsintime(ii,2)) < 0.0001)));
end



% AdWell fine mesh
filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive\181005_4';

fullpath = [filepath_simdata '\solution\monitors\fluid-superficial-x-velocity.txt'];
[time, U_f] = importMonitoredQuantities(fullpath);

fullpath = [filepath_simdata '\solution\monitors\solid-superficial-x-velocity.txt'];
[time, U_s] = importMonitoredQuantities(fullpath);

plot(time,U_f,'-','Color',[0 112/256 192/256],'LineWidth',1);
plot(time,U_s,'-','Color',[204/256 102/256 0],'LineWidth',1);


for ii = 1:size(pointsintime,1)
    mean_U_f_fine(ii) = mean(U_f(find(abs(time-pointsintime(ii,1)) < 0.0001):find(abs(time-pointsintime(ii,2)) < 0.0001)));
    mean_U_s_fine(ii) = mean(U_s(find(abs(time-pointsintime(ii,1)) < 0.0001):find(abs(time-pointsintime(ii,2)) < 0.0001)));
end

TightFigure(2);


for ii = 1:size(pointsintime,1)
    ratio = mean_U_f_coarse(ii)/mean_U_f_fine(ii);
    text(mean(pointsintime(ii,:)),mean_U_f_fine(ii)+0.08,['r_f = ' num2str(ratio,3)],'Color','black','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom');
    ratio = mean_U_s_coarse(ii)/mean_U_s_fine(ii);
    text(mean(pointsintime(ii,:)),mean_U_s_fine(ii)+0.12,['r_s = ' num2str(ratio,3)],'Color','black','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom');
end

set(gca,...
    'FontSize',12,...
    'xLim',[0 16],...
    'yLim',[0 1.8]);

arrowheight = 1.65
arrow([12.6 arrowheight], [11 arrowheight])
arrow([14.4 arrowheight], [16 arrowheight])

arrowheight = 0.55
arrow([12.6 arrowheight], [11 arrowheight])
arrow([14.4 arrowheight], [16 arrowheight])


legend('U_f (Intermediate)', 'U_s (Intermediate)', 'U_f (Fine)', 'U_s (Fine)', 'Location', 'Northwest')

% Print to files
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Wellbore\Figures\';
fig_name = 'mesh_mp';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 

