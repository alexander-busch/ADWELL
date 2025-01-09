function [ subfig ] = CreateSubplot( subplot_definition, subplot_title, label_x, label_y, sub_fontsize )
%	CreateSubplot Format subplot
%   Detailed explanation goes here

subfig = subplot(subplot_definition(1), subplot_definition(2), subplot_definition(3));
hold on;
grid on;
box on;

title(subplot_title,'FontSize',sub_fontsize);
xlabel(label_x,'FontSize',sub_fontsize);
ylabel(label_y,'FontSize',sub_fontsize);

set(gca,...
    'FontSize',sub_fontsize)

% set(gca,...
%     'XScale','log',...
%     'YScale','log');

end
