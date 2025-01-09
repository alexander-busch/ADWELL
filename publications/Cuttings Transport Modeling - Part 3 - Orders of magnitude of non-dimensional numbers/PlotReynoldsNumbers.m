% Plot Reynolds Numbers

%% Reynolds number for Ostwald fluid

fig_title = 'Generalized Reynolds number for Ostwald (PL) behavior';

CreateFigure( fig_title,...
    'Hydraulic diameter d_h [m]',...
    'Generalized Reynolds number \rho_f U d_h / \eta_f [-]',...
    {'lin' 'log'}, 'PowerPoint');
set(gca,'ColorOrder', ColorSet);
% Loop all fluids
plot_h = zeros(length(PV)+3,1);
legendentries = cell(length(PV)+3,1);
% ii=1; ii=2; ii=3; ii=4; ii=5;
for ii=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % jj=1; jj=2;    
    for jj=1:max(size(d_h))
        % Pipe
        plot_h(ii) = scatter([d_h(jj,1) d_h(jj,1)],Re_PL(jj,1:2,ii),'o','filled');
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,1) d_h(jj,1)],Re_PL(jj,1:2,ii));
        set(gca,'ColorOrderIndex',COI);
        % Annulus
        scatter([d_h(jj,2) d_h(jj,2)],[Re_PL(jj,3,ii) Re_PL(jj,6,ii)],'o','filled');
        set(gca,'ColorOrderIndex',COI);
        scatter([d_h(jj,2) d_h(jj,2)],[Re_PL(jj,3,ii) Re_PL(jj,6,ii)],'wo','MarkerFaceColor', 'w','SizeData', 5);
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,2) d_h(jj,2)],[Re_PL(jj,3,ii) Re_PL(jj,6,ii)]);
        set(gca,'ColorOrderIndex',COI);
    end
    % Critical Re
    plot(d_h(:,2),ones(length(d_h),1)*Re_cr_PL(ii,1),'-');
    set(gca,'ColorOrderIndex',COI);
%     plot(d_h(:,2),ones(length(d_h),1)*Re_cr_PL(ii,2),'--');
%     set(gca,'ColorOrderIndex',COI);
    % Legend (Fluid colors)
    plot_h(ii) = scatter(100*[d_h(jj,1) d_h(jj,1)],Re_PL(jj,1:2,ii),'s','filled');
    legendentries{ii} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),1), ', PV = ', num2str(PV(ii)),')'];
end
% Legend
plot_h(ii+1) = scatter(10,10,'ko','filled'); % Pipe
legendentries(ii+1)={'Pipe'};
plot_h(ii+2) = scatter(10,10,'ko'); % Annulus
legendentries(ii+2)={'Annulus'};
plot_h(ii+3) = plot(10,10,'k-'); % Re_cr_PL = f(n_P)
legendentries(ii+3)={'Re_c_r = f(n_P)'};
% plot_h(ii+4) = plot(10,10,'k--'); % Re_cr_PL = f(n_A)
% legendentries(ii+4)={'Re_c_r = f(n_A)'};
legend(plot_h,legendentries,'Location','northwest');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'xlim',[0 1],...
    'ylim',[1e2 2e6]);
TightFigure(gca);
legend(plot_h,legendentries,'Location','northwest');legend off

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'




CreateFigure( [fig_title ', fluid ' num2str(fluid)],...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]',...
    'Generalized Reynolds number \rho_f U d_h / \eta_f [-]',...
    {'lin' 'log'}, 'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[Re_PL(:,2,fluid) Re_PL(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[Re_PL(:,1,fluid) Re_PL(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_PL(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_PL(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_PL(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_PL(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot Re_cr
plot(sections,ones(length(sections),1)*Re_cr_PL(fluid,1),'k-');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e2 1e5]);
legend('Drll pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Re_c_r = f (n_P)',...
    'Location','northeast');
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates'...
    ', recommended drill pipe rotation rates',...
    ' and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2),...
    ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,1.15*max(ylim),textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Reynolds number for Bingham fluid

fig_title = 'Generalized Reynolds number for Bingham (YP/PV) behavior';

CreateFigure( fig_title,...
    'Hydraulic diameter d_h [m]',...
    'Generalized Reynolds number \rho_f U d_h / \eta_f [-]',...
    {'lin' 'log'}, 'PowerPoint');
set(gca,'ColorOrder', ColorSet);
% Loop all fluids
plot_h = zeros(length(PV)+4,1);
legendentries = cell(length(PV)+4,1);
% ii=1; ii=2; ii=3; ii=4; ii=5;
for ii=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % jj=1; jj=2;
    for jj=1:max(size(d_h))
        % Pipe
        plot_h(ii) = scatter([d_h(jj,1) d_h(jj,1)],Re_YPPV(jj,1:2,ii),'o','filled');
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,1) d_h(jj,1)],Re_YPPV(jj,1:2,ii));
        set(gca,'ColorOrderIndex',COI);
        % Annulus
        scatter([d_h(jj,2) d_h(jj,2)],[Re_YPPV(jj,3,ii) Re_YPPV(jj,6,ii)],'o','filled');
        set(gca,'ColorOrderIndex',COI);
        scatter([d_h(jj,2) d_h(jj,2)],[Re_YPPV(jj,3,ii) Re_YPPV(jj,6,ii)],'wo','MarkerFaceColor', 'w','SizeData', 5);
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,2) d_h(jj,2)],[Re_YPPV(jj,3,ii) Re_YPPV(jj,6,ii)]);
        set(gca,'ColorOrderIndex',COI);
    end
    % Critical Re    
    plot(d_h(:,1),ones(length(d_h),1)*Re_cr_PL(ii,1),'-');
    set(gca,'ColorOrderIndex',COI);
    plot(d_h(:,2),ones(length(d_h),1)*Re_cr_PL(ii,2),'--');
    set(gca,'ColorOrderIndex',COI);
    % Legend (Fluid colors)
    plot_h(ii) = scatter(100*[d_h(jj,1) d_h(jj,1)],Re_PL(jj,1:2,ii),'s','filled');
    legendentries{ii} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),1), ', PV = ', num2str(PV(ii)),')'];
end
% Legend
plot_h(ii+1) = scatter(10,10,'ko','filled'); % Pipe
legendentries(ii+1)={'Pipe'};
plot_h(ii+2) = scatter(10,10,'ko'); % Annulus
legendentries(ii+2)={'Annulus'};
plot_h(ii+3) = plot(10,10,'k-'); % Re_cr_PL = f(n_P)
legendentries(ii+3)={'Re_c_r = f(n_P)'};
plot_h(ii+4) = plot(10,10,'k--'); % Re_cr_PL = f(n_A)
legendentries(ii+4)={'Re_c_r = f(n_A)'};
% legend(plot_h,legendentries,'Location','northwest');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'xlim',[0 1],...
    'ylim',[1e2 2e6]);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



CreateFigure( [fig_title ', fluid ' num2str(fluid)],...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]',...
    'Generalized Reynolds number \rho_f U d_h / \eta_f [-]',...
    {'lin' 'log'}, 'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[Re_YPPV(:,2,fluid) Re_YPPV(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[Re_YPPV(:,1,fluid) Re_YPPV(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_YPPV(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_YPPV(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_YPPV(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_YPPV(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot Re_cr
plot(sections,ones(length(sections),1)*Re_cr_PL(fluid,1),'-','Color',col_dp);
plot(sections,ones(length(sections),1)*Re_cr_PL(fluid,2),'--','Color',col_an1);

legend('Drll pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Re_c_r = f (n_P)',...
	'Re_c_r = f (n_A)',...
    'Location','northeast');
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e2 1e5]);
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates'...
    ', recommended drill pipe rotation rates',...
    ' and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2),...
    ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,1.15*max(ylim),textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'




%% Classical Reynolds number based on Bingham viscosity and Hedström number (after Wilson and Thomas (2006))
fig_title = 'Classic Bingham Reynolds number (PV based) for Bingham (YP/PV) behavior';

CreateFigure( fig_title,...
    'Hydraulic diameter d_h [m]',...
    'Bingham Reynolds number \rho_f U d_h / \mu_P_V[-]',...
    {'lin' 'log'}, 'PowerPoint');
set(gca,'ColorOrder', ColorSet);
% Loop all fluids
% ii=1; ii=2; ii=3; ii=4; ii=5;
plot_h = zeros(length(PV)+4,1);
legendentries = cell(length(PV)+4,1);
for ii=1:length(PV)
    % Color index
    COI = get(gca,'ColorOrderIndex');
    % Loop d_h and plot pointwise in order to express range with endpoints
    % and a line in between
    % jj=1; jj=2;
    for jj=1:max(size(d_h))
                % Pipe
        plot_h(ii) = scatter([d_h(jj,1) d_h(jj,1)],Re_Bi(jj,1:2,ii),'o','filled');
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,1) d_h(jj,1)],Re_Bi(jj,1:2,ii));
        set(gca,'ColorOrderIndex',COI);
        % Annulus
        scatter([d_h(jj,2) d_h(jj,2)],[Re_Bi(jj,3,ii) Re_Bi(jj,6,ii)],'o','filled');
        set(gca,'ColorOrderIndex',COI);
        scatter([d_h(jj,2) d_h(jj,2)],[Re_Bi(jj,3,ii) Re_Bi(jj,6,ii)],'wo','MarkerFaceColor', 'w','SizeData', 5);
        set(gca,'ColorOrderIndex',COI);
        plot([d_h(jj,2) d_h(jj,2)],[Re_Bi(jj,3,ii) Re_Bi(jj,6,ii)]);
        set(gca,'ColorOrderIndex',COI);
    end
    % Critical Re    
    plot(d_h(:,1),Re_cr_Bi(:,1,ii),'-');
    set(gca,'ColorOrderIndex',COI);
    plot(d_h(:,2),Re_cr_Bi(:,3,ii),'--');
    set(gca,'ColorOrderIndex',COI);
    % Legend (ii colors)
    plot_h(ii) = scatter(100*[d_h(jj,1) d_h(jj,1)],Re_PL(jj,1:2,ii),'s','filled');
    legendentries{ii} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),1), ', PV = ', num2str(PV(ii)),')'];
    % Add information on fluids and critical Re
    if ii<4
        string = '';
        stringcolor = 'k';
    elseif ii==4
        string = 'Fluid 1-4';
        stringcolor = 'k';
        plot_h(length(plot_h)) = plot(d_h(:,2),Re_cr_Bi(:,3,ii),'k--');
        legendentries(length(legendentries))={'Re_c_r = f(He)'};
    else
        string = ['Fluid ' num2str(ii)];
        stringcolor = ColorSet(COI,:);
    end
    text(max(d_h(:,2))+0.01,max(Re_cr_Bi(:,3,ii)),string,...
        'Color',stringcolor,...
        'FontSize',10,...
        'FontWeight','normal',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle',...
        'Rotation',0);
    if ii==8
        string = 'PV = 15, Fluids 2, 5, 8';
    elseif ii==9
        string = 'PV = 30, Fluids 4, 6, 9';
    elseif ii==10
        string = 'PV = 50, Fluids 5, 7, 10';
    else
        string = '';
    end
    text(0.5,Re_Bi(3,6,ii),string,...
        'Color',ColorSet(COI,:),...
        'FontSize',10,...
        'FontWeight','normal',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle',...
        'Rotation',0);
end
% Legend
plot_h(ii+1) = scatter(10,10,'ko','filled'); % Pipe
legendentries(ii+1)={'Pipe'};
plot_h(ii+2) = scatter(10,10,'ko'); % Annulus
legendentries(ii+2)={'Annulus'};
plot_h(ii+3) = plot(10,10,'k-'); % Re_cr_Bi = f(n_P)
legendentries(ii+3)={'Re_c_r = f(He) (Pipe)'};
plot_h(ii+4) = plot(10,10,'k--'); % Re_cr_Bi = f(n_A)
legendentries(ii+4)={'Re_c_r = f(He) (Annulus)'};
% legend(plot_h,legendentries,'Location','north');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'xlim',[0 1],...
    'ylim',[1e2 2e6]);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



CreateFigure( [fig_title ', fluid ' num2str(fluid)],...
    'Hydraulic diameter [m], Drill pipe diameter [in], Wellbore section [in]',...
    'Bingham Reynolds number \rho_f U d_h / \mu_P_V[-]',...
    {'lin' 'log'}, 'PowerPoint');
% Plot zero rotation in background
bar_handle = bar(sections,[Re_Bi(:,2,fluid) Re_Bi(:,6,fluid)],1.0);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[Re_Bi(:,1,fluid) Re_Bi(:,3,fluid)],1.0);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot rec rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_Bi(:,7,fluid)],0.66);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an2; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_Bi(:,4,fluid)],0.66);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an2;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot max rpm in foreground
bar_handle = bar(sections,[zeros(length(d_h),1) Re_Bi(:,8,fluid)],0.33);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an3; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
bar_handle = bar(sections,[zeros(length(d_h),1) Re_Bi(:,5,fluid)],0.33);
bar_handle(1).FaceColor = 'w'; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an3;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% Plot Re_cr
plot(sections,Re_cr_Bi(:,1,fluid),'-','Color',col_dp);
plot(sections,Re_cr_Bi(:,3,fluid),'--','Color',col_an1);
% Legend
legend('Drll pipe',...
    'Annulus, 0 rpm',...
    'Annulus, 1 rpm',...
    'Annulus, 2 rpm',...
    'Re_c_r = f(He) (Pipe)',...
    'Re_c_r = f(He) (Annulus)',...
    'Location','southeast');
% Format
set(gca,...
    'Xdir', 'reverse',...
    'ylim',[1e2 1e5]);
xticklabels(flip(matrix));
textstring = ['For respective min. and max. bulk flow rates'...
    ', recommended drill pipe rotation rates',...
    ' and fluid ' num2str(fluid) ' (n = ' num2str(n_A(fluid),2),...
    ', K = ' num2str(n_A(fluid),2) ').'];
% text(length(d_h)+0.5,1.15*max(ylim),textstring);
TightFigure(gca);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
