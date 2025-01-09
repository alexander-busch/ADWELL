% Plot apparent viscosities at the wall vs. hole sections


%% Apparent viscosity based on true wall shear rate for Ostwald fluid

CreateFigure( 'Apparent viscosity for Ostwald fluid (PL) behavior',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Apparent viscosity [Pa.s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[eta_f_PL(:,1,fluid) eta_f_PL(:,3,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[eta_f_PL(:,2,fluid) eta_f_PL(:,6,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_PL(:,4,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_PL(:,7,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_PL(:,5,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_PL(:,8,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e-2 1e1]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Location','northeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates, recommended drill pipe rotation rates and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2) ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,max(ylim)+2,textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Apparent viscosity based on true wall shear rate for Bingham fluid

CreateFigure( 'Apparent viscosity for Bingham fluid (YP/PV) behavior',...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]', 'Apparent viscosity [Pa.s]',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[eta_f_YPPV(:,1,fluid) eta_f_YPPV(:,3,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[eta_f_YPPV(:,2,fluid) eta_f_YPPV(:,6,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_YPPV(:,4,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_YPPV(:,7,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_YPPV(:,5,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) eta_f_YPPV(:,8,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e-2 1e1]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Location','northeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates, recommended drill pipe rotation rates and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2) ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,max(ylim)+2,textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
