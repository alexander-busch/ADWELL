%% Shields number for Ostwald fluid
fig_title = 'Shields number for Ostwald (PL) behavior';

CreateFigure( [fig_title ', and d_p = ' num2str(d_p(particle)*1000,2) ' mm'],...
    'Hydraulic diameter d_h [m]',...
    'Normalized Shields number Sh/Sh_c_r',...
    {'lin' 'log'},'PowerPoint');
set(gca,'ColorOrder', ColorSet);
% Loop all fluids
plot_h = zeros(length(PV)+2,1);
legendentries = cell(length(PV)+2,1);
% ii=1; ii=2; ii=3; ii=4; ii=5;
for ii=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % jj=1; jj=2;    
    for jj=1:max(size(d_h))
        % Pipe
        scatter([d_h(jj,1) d_h(jj,1)],Sh_PL(jj,1:2,ii),'o','filled');
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,1) d_h(jj,1)],Sh_PL(jj,1:2,ii));
        set(gca,'ColorOrderIndex',COI);
        % Annulus
        scatter([d_h(jj,2) d_h(jj,2)],[Sh_PL(jj,3,ii) Sh_PL(jj,6,ii)],'o','filled');
        set(gca,'ColorOrderIndex',COI);
        scatter([d_h(jj,2) d_h(jj,2)],[Sh_PL(jj,3,ii) Sh_PL(jj,6,ii)],'wo','MarkerFaceColor', 'w','SizeData', 5);
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,2) d_h(jj,2)],[Sh_PL(jj,3,ii) Sh_PL(jj,6,ii)]);
        set(gca,'ColorOrderIndex',COI);
    end
    % Critical Sh
    plot(d_h(:,2),ones(length(d_h),1),'k-');
    set(gca,'ColorOrderIndex',COI);
    % Legend (Fluid colors)
    plot_h(ii) = scatter(100*[d_h(jj,1) d_h(jj,1)],Sh_PL(jj,1:2,ii),'s','filled');
    legendentries{ii} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ', PV = ', num2str(PV(ii),2),')'];
end
% Legend
plot_h(ii+1) = scatter(10,10,'ko'); % Annulus
legendentries(ii+1)={'Annulus'};
plot_h(ii+2) = plot(d_h(:,2),ones(length(d_h),1),'k-');
legendentries(ii+2)={'Sh_c_r'};
% legend(plot_h,legendentries,'Location','southwest');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'xlim',[0 1],...
    'ylim',[1e-3 1e1]);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'

CreateFigure( [fig_title ', fluid ' num2str(fluid) ', and d_p = ' num2str(d_p(particle)*1000,2) ' mm'],...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]',...
    'Normalized Shields number Sh/Sh_c_r',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[Sh_PL(:,2,fluid) Sh_PL(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[Sh_PL(:,1,fluid) Sh_PL(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_PL(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_PL(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_PL(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_PL(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot Sh_cr
plot(sections,ones(length(sections),1),'k-');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e-1 1e1]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Sh_c_r',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates'...
    ', recommended drill pipe rotation rates',...
    ' fluid ' num2str(fluid) ' (YP = ' num2str(YP(fluid),2),...
    ', PV = ' num2str(PV(fluid),2),...
    ' and d_p = ' num2str(d_p(particle)*1000,2) ' mm).'];
% text(length(d_h)+0.5,1.15*max(ylim),textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Shields number for Bingham fluid

CreateFigure( [fig_title ', and d_p = ' num2str(d_p(particle)*1000,2) ' mm'],...
    'Hydraulic diameter d_h [m]',...
    'Normalized Shields number Sh/Sh_c_r',...
    {'lin' 'log'},'PowerPoint');
set(gca,'ColorOrder', ColorSet);
% Loop all fluids
plot_h = zeros(length(PV)+2,1);
legendentries = cell(length(PV)+2,1);
% ii=1; ii=2; ii=3; ii=4; ii=5;
for ii=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % jj=1; jj=2;    
    for jj=1:max(size(d_h))
        % Pipe
        scatter([d_h(jj,1) d_h(jj,1)],Sh_YPPV(jj,1:2,ii),'o','filled');
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,1) d_h(jj,1)],Sh_YPPV(jj,1:2,ii));
        set(gca,'ColorOrderIndex',COI);
        % Annulus
        scatter([d_h(jj,2) d_h(jj,2)],[Sh_YPPV(jj,3,ii) Sh_YPPV(jj,6,ii)],'o','filled');
        set(gca,'ColorOrderIndex',COI);
        scatter([d_h(jj,2) d_h(jj,2)],[Sh_YPPV(jj,3,ii) Sh_YPPV(jj,6,ii)],'wo','MarkerFaceColor', 'w','SizeData', 5);
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,2) d_h(jj,2)],[Sh_YPPV(jj,3,ii) Sh_YPPV(jj,6,ii)]);
        set(gca,'ColorOrderIndex',COI);
    end
    % Critical Sh
    plot(d_h(:,2),ones(length(d_h),1),'k-');
    set(gca,'ColorOrderIndex',COI);
    % Legend (Fluid colors)
    plot_h(ii) = scatter(100*[d_h(jj,1) d_h(jj,1)],Sh_YPPV(jj,1:2,ii),'s','filled');
    legendentries{ii} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ', PV = ', num2str(PV(ii),2),')'];
end
% Legend
plot_h(ii+1) = scatter(10,10,'ko'); % Annulus
legendentries(ii+1)={'Annulus'};
plot_h(ii+2) = plot(d_h(:,2),ones(length(d_h),1),'k-');
legendentries(ii+2)={'Sh_c_r'};
% legend(plot_h,legendentries,'Location','southwest');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'xlim',[0 1],...
    'ylim',[1e-3 1e1]);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


fig_title = 'Shields number for Bingham (YP/PV) behavior';

CreateFigure( [fig_title ', fluid ' num2str(fluid) ', and d_p = ' num2str(d_p(particle)*1000,2) ' mm'],...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]',...
    'Normalized Shields number Sh/Sh_c_r',...
    {'lin' 'log'},'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[Sh_YPPV(:,2,fluid) Sh_YPPV(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[Sh_YPPV(:,1,fluid) Sh_YPPV(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_YPPV(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_YPPV(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_YPPV(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Sh_YPPV(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot Sh_cr
plot(sections,ones(length(sections),1),'k-');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e-1 1e1]);
legend('Drill pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Sh_c_r',...
    'Location','southeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates'...
    ', recommended drill pipe rotation rates',...
    ' fluid ' num2str(fluid) ' (YP = ' num2str(YP(fluid),2),...
    ', PV = ' num2str(PV(fluid),2),...
    ' and d_p = ' num2str(d_p(particle)*1000,2) ' mm).'];
% text(length(d_h)+0.5,1.15*max(ylim),textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'

