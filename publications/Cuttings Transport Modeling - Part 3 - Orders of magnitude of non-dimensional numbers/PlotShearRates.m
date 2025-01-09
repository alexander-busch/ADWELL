% Plot Newtonian and PL shear rates vs. hole sections


%% Newtonian wall shear rate
CreateFigure( 'Newtonian wall shear rate for annulus and drill pipe as well as recommended drill pipe rotation rates vs. hole section',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Newtonian wall shear rate [1/s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[SR_N(:,2) SR_N(:,6)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[SR_N(:,1) SR_N(:,3)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_N(:,7)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_N(:,4)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_N(:,8)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_N(:,5)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1 1400]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = 'For respective min. and max. bulk flow rates and recommended drill pipe rotation rates.';
% text(l'ength(d_h)+0.5,max(ylim)+200,textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Limiting PL wall shear rate
CreateFigure( 'Limiting PL wall shear rate for annulus and drill pipe as well as recommended drill pipe rotation rates vs. hole section',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Wall shear rate [1/s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[SR_PL(:,2,fluid) SR_PL(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[SR_PL(:,1,fluid) SR_PL(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_PL(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_PL(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_PL(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_PL(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1 1400]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates, recommended drill pipe rotation rates and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2) ', K = ' num2str(n_A(fluid),2) ').'];
% text(l'ength(d_h)+0.5,max(ylim)+200,textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'

%% Limiting YP/PV wall shear rate
CreateFigure( 'Limiting YP/PV wall shear rate for annulus and drill pipe as well as recommended drill pipe rotation rates vs. hole section',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Wall shear rate [1/s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[SR_YPPV(:,2,fluid) SR_YPPV(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[SR_YPPV(:,1,fluid) SR_YPPV(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_YPPV(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_YPPV(:,4,fluid)],0.66);                                                 
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) SR_YPPV(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) SR_YPPV(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1 1400]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates, recommended drill pipe rotation rates and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2) ', K = ' num2str(n_A(fluid),2) ').'];
% text(l'ength(d_h)+0.5,max(ylim)+200,textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
