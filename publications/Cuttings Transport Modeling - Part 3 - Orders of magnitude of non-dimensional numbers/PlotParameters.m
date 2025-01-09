% Plot parameters vs. hole sections


%% Geometry vs. hole section
CreateFigure( 'Hole, drill pipe and hydraulic diameter vs. hole section',...
    'Wellbore section [in]', 'Hole, inner and outer drill pipe, and hydraulic diameter [in]',...
    {'lin' 'lin'},'PowerPoint');
bar_handle = bar(sections,[d_o(:,1)/in2m d_o(:,2)/in2m],0.9);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = col_dp;
bar_handle(2).FaceColor = col_an1; bar_handle(2).EdgeColor = col_an1;
bar_handle = bar(sections,[d_i(:,1)/in2m d_i(:,2)/in2m],0.9);
bar_handle(1).FaceColor = col_dp; bar_handle(1).EdgeColor = 'k';
bar_handle(2).FaceColor = 'w'; bar_handle(2).EdgeColor = col_an1;
set(get(get(bar_handle(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'Xdir', 'reverse');
legend('Drill pipe','Annulus');
xticklabels(flip(matrix2(:,3)));
TightFigure(gca); % Handles to axis, if only one just provide one
for ii=1:length(d_h)
    text(ii+0.15,d_o(length(d_h)+1-ii,2)/in2m,...
        num2str(d_o(length(d_h)+1-ii,2)/in2m,2),'Color',col_an1,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    text(ii-0.15,d_o(length(d_h)+1-ii,1)/in2m,...
        num2str(d_o(length(d_h)+1-ii,1)/in2m,2),'Color',col_dp,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    text(ii+0.15,d_i(length(d_h)+1-ii,2)/in2m,...
        num2str(d_i(length(d_h)+1-ii,2)/in2m,2),'Color',...
        col_an1,'HorizontalAlignment','Center','VerticalAlignment','Top');
    text(ii+0.15,d_i(length(d_h)+1-ii,2)/in2m + d_h(length(d_h)+1-ii,2)/2/in2m,...
        num2str(d_h(length(d_h)+1-ii,2)/in2m,2),'Color',...
        'w','HorizontalAlignment','Center','VerticalAlignment','Middle');
    text(ii-0.15,d_h(length(d_h)+1-ii,1)/2/in2m,...
        num2str(d_h(length(d_h)+1-ii,1)/in2m,1),'Color',...
        'w','HorizontalAlignment','Center','VerticalAlignment','Middle');
end

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'




%% Fluid flow rate vs. hole section
CreateFigure( 'Drilling fluid flow rate vs. hole section',...
    'Wellbore section [in]', 'Drilling fluid flow rate [L/min]',...
    {'lin' 'lin'},'PowerPoint');
bar_handle(1) = bar(sections,Q(:,2)/gpm2m3ps*gpm2Lpm,0.8,'FaceColor',lightblue,'EdgeColor',lightblue);
bar_handle(2) = bar(sections,Q(:,1)/gpm2m3ps*gpm2Lpm,0.8,'FaceColor','w','EdgeColor',lightblue);
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'Xdir', 'reverse');
legend('Recommended flow rate (K&M)');
xticklabels(flip(matrix2(:,3)));
% Minimize white space of figure 
TightFigure(gca); % Handles to axis, if only one just provide one
for ii=1:length(d_h)
    text(ii,Q(length(d_h)+1-ii,2)/gpm2m3ps*gpm2Lpm,...
        num2str(Q(length(d_h)+1-ii,2)/gpm2m3ps*gpm2Lpm,3),'Color',lightblue,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    text(ii,Q(length(d_h)+1-ii,1)/gpm2m3ps*gpm2Lpm,...
        num2str(Q(length(d_h)+1-ii,1)/gpm2m3ps*gpm2Lpm,3),'Color',lightblue,...
        'HorizontalAlignment','Center','VerticalAlignment','Top');
end

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Drill pipe rotation vs. hole section
CreateFigure( 'Drill pipe rotation rate vs. hole sections',...
    'Wellbore section [in]', 'Drill pipe rotation rate [rpm]',...
    {'lin' 'lin'},'PowerPoint');
set(gca,'ylim',[0 500]);
bar_handle(1) = bar(sections,rpm(:,2),0.8,'FaceColor',purple,'EdgeColor',purple);
bar_handle(2) = bar(sections,rpm(:,1),0.8,'FaceColor','w','EdgeColor',purple);
set(get(get(bar_handle(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'Xdir', 'reverse');
legend('Recommended rotation rate (K&M)');
xticklabels(flip(matrix2(:,3)));
% Minimize white space of figure 
TightFigure(gca); % Handles to axis, if only one just provide one
for ii=1:length(d_h)
    text(ii,rpm(length(d_h)+1-ii,2),...
        num2str(rpm(length(d_h)+1-ii,2),3),'Color',purple,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    text(ii,rpm(length(d_h)+1-ii,1),...
        num2str(rpm(length(d_h)+1-ii,1),3),'Color',purple,...
        'HorizontalAlignment','Center','VerticalAlignment','Top');
end


% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
