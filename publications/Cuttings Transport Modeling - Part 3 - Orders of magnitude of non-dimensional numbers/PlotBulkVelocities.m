% Plot bulk velocities


CreateFigure( 'Bulk fluid velocity magnitudes for annulus and drill pipe as well as recommended drill pipe rotation rates vs.  hole section',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Bulk fluid velocity magnitude [m/s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rpm in background
bar_handle = bar(sections,[U(:,2) U(:,6)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[U(:,1) U(:,3)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot recommended rpm in center
bar_handle = bar(sections,[zeros(length(d_h),1) U(:,7)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) U(:,4)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) U(:,8)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) U(:,5)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'YScale','log',...
    'ylim',[0.1 10]);
legend('Drill pipe',...
    'Annulus (0 x rpm_r)',...
    'Annulus (1 x rpm_r)',...
    'Annulus (2 x rpm_r)',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = 'For respective min. and max. bulk flow rates and recommended drill pipe rotation rates.';
% text(length(d_h)+0.5,max(ylim)+1,textstring);
TightFigure(gca);


% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'