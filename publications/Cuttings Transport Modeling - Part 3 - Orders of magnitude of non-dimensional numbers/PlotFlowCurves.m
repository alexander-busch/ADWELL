% Plot flow curves

% Plot PAC? 
PAC = 'Y'; % 'N'
if strcmp(PAC,'Y')
    % Cross and Carreau fits of rheometric data
    p = pwd;
    load([p(1) ':\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\AdWell\Particle settling\Particle settling in a shear flow\Cross.mat']);
    load([p(1) ':\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\AdWell\Particle settling\Particle settling in a shear flow\Carreau.mat']);
    clear p;
    % Fluids
    fluidlist = {'H2O' ...
        'PAC1' ...
        'PAC2' ...
        'PAC3' ...
        'PAC4'};
end

% fig = CreateFigure( '', '', '', {'log','log'}, 'PowerPoint' );
% axis off;
% fig_sub = tight_subplot(1,2,[.1 .06],[.06 .03],[.06 .02]);
% 
% subplot(222); set(gca,'Visible','Off');

%% Shear stress vs. shear rate, PL fit of lower shear rate range

PH_FC(1) = CreateFigure(...
    'Shear stress as a function shear rate w/ PL fit based on typical annular shear rate range',...
    'Shear rate [1/s]', 'Shear stress [Pa]', {'log' 'log'},'PowerPoint');
set(gca,'ColorOrder', ColorSet,'xlim',[1e0 1200], 'ylim',[1e-3 1e2]);

% subplot(232); hold on;
% set(gca,'ColorOrder', ColorSet,...
%     'XScale','log','YScale','log',...
%     'xlim',[1e0 1200],'ylim',[1e-3 1e2]);

% Shaded shear rate intervalls for pipe and annular flow
SRI = [area([SR_Fann(4) SR_Fann(6)],[max(ylim) max(ylim)])...
area([SR_Fann(2) SR_Fann(1)],[max(ylim) max(ylim)])];
    for aa = 1:length(SRI)
        SRI(aa).FaceColor = [0.8 0.8 0.8];
        SRI(aa).EdgeColor = 'none';
        SRI(aa).FaceAlpha = 0.5;
    end
% Initialize arrays
tau_YPPV=zeros(length(YP),length(SR_range));
tau_PL_A=zeros(length(YP),length(SR_range));
tau_PL_P=zeros(length(YP),length(SR_range));
legendentries=cell(length(YP),2);
fluidtable=cell(length(YP),7);
plot_h=zeros(length(YP),2);
% Plot shear stress of representative drilling fluids
for ii=1:length(YP)
    COI = get(gca,'ColorOrderIndex');
    % YP/PV
    tau_YPPV(ii,:) = (YP(ii)+PV(ii).*SR_range)';
    plot_h(ii,1)=plot(SR_range,tau_YPPV(ii,:),':');
    % PL
    tau_PL_A(ii,:) = (K_A(ii).*SR_range.^(n_A(ii)))';
    set(gca,'ColorOrderIndex',COI);
    plot_h(ii,2)=plot(SR_range,tau_PL_A(ii,:),'--');
    % Reference values at 3 and 100 Fann rpm
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(6),tau_Fann(ii,6),'o');
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(4),tau_Fann(ii,4),'o');
    legendentries{ii,1} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),')'];
    legendentries{ii,2} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),...
        ' \rightarrow n = ', num2str(n_A(ii),2), ', K = ', num2str(K_A(ii),3),')'];
    fluidtable{ii,1} = num2str(ii);
    fluidtable{ii,2} = num2str(YP(ii),2);
    fluidtable{ii,3} = num2str(PV(ii),2);
    % Pipe to be filled later
    fluidtable{ii,6} = num2str(n_A(ii),2);
    fluidtable{ii,7} = num2str(K_A(ii),2);

    
end
% Plot PAC
if strcmp(PAC,'Y')
    for ii=1:length(fluidlist)
        % Plot Cross fit
        eta = Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*SR_range).^Cross.n{ii});       
        tau = eta.*SR_range;
        plot(SR_range,tau,'-','Color','k');
        % Plot Carreau fit
%         plot(SR,Carreau.mu_inf{ii}+(Carreau.mu_0{ii}-Carreau.mu_inf{ii}).*((1+(Carreau.lambda{ii}.*SR_range).^2).^((Carreau.n{ii}-1)./2)),'-','Color','k')
        pos = min(xlim)+10^(log10(min(xlim))-1);
        text(pos,....
            Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*pos).^Cross.n{ii})*pos,fluidlist{ii},...
            'HorizontalAlignment','Left',...
            'VerticalAlignment','Bottom',...
            'Rotation',25);
    end    
end
plot_h = [plot_h(:,1); plot_h(:,2)];
legendentries = [legendentries(:,1); legendentries(:,2)];
% legend_h = legend(plot_h,legendentries,'Location','southeast','FontSize',10); %eastoutside
% % legend_h = legend(plot_h(:,2),legendentries(:,3),'Location','southeast','FontSize',10);
% Plot text
% text(0.6,1,'YP [lbf/100ft²]\newlinePV [mPa.s]\newlinen [-]\newlineK[mPa.s^n]');
% colors = get(groot,'defaultAxesColorOrder');
tx = 30;
ty = 1.5e-3;
text(tx,ty,...
    {'Typical','annular flow','shear rate range'},...
    'Color','w',...
    'FontSize',18,...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'Rotation',0);
% Second set of axis
ax1 = gca; % current axes
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top','YAxisLocation','right',...   
    'Xlim',xlim/1.703,'Ylim',ylim/0.4788026,...
    'XScale','log','YScale','log',...
    'Xtick',[3 6 100 200 300 600],...
    'Color','none');
xlabel('Fann viscometer speed [rpm]'); ylabel('Shear stress [lbf/100ft²]');
% Minimize white space of figure 
TightFigure([ax1 ax2]); % Handles to axis, if only one just provide one

% 
% % Get positions
% pos_ax1=get(ax1,'position'); 
% pos_ax2=get(ax2,'position'); 
% pos_fig=get(gcf,'position');
% % % pos_leg=get(legend_h,'position');
% 
% % Resize factor
% factor = 1+% % pos_leg(3);
% 
% % Resize axis
% pos_ax1(1)=pos_ax1(1)/factor; pos_ax1(3)=pos_ax1(3)/factor-0.015;
% pos_ax2(1)=pos_ax2(1)/factor; pos_ax2(3)=pos_ax2(3)/factor-0.015;
% 
% % Resize figure
% pos_fig(3)=factor*pos_fig(3)+0.015;
% set(gcf,'position',pos_fig);
% 
% % Reset axis
% set(ax1,'position',pos_ax1);
% set(ax2,'position',pos_ax2);
% 
% % Programatically move the Legend
% % % pos_leg(3)=% % pos_leg(3)/factor; % % pos_leg(1)=1-% % pos_leg(3)-0.015;
% % set(legend_h,'Position',% % pos_leg);


% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'




%% Shear stress vs. shear rate, PL fit of higher shear rate range
PH_FC(2) = CreateFigure(...
    'Shear stress as a function shear rate w/ PL fit based on typical pipe shear rate range',...
    'Shear rate [1/s]', 'Shear stress [Pa]', {'lin' 'lin'},'PowerPoint');
set(gca,'ColorOrder', ColorSet,'xlim',[0 1200], 'ylim',[0 80]);
% Shaded shear rate intervalls for pipe and annular flow
SRI = [area([SR_Fann(4) SR_Fann(6)],[max(ylim) max(ylim)])...
area([SR_Fann(2) SR_Fann(1)],[max(ylim) max(ylim)])];
    for aa = 1:length(SRI)
        SRI(aa).FaceColor = [0.8 0.8 0.8];
        SRI(aa).EdgeColor = 'none';
        SRI(aa).FaceAlpha = 0.5;
    end

legendentries=cell(length(YP),2);
plot_h=zeros(length(YP),2);
% Plot shear stress of representative drilling fluids
for ii=1:length(YP)
    COI = get(gca,'ColorOrderIndex');
    % YP/PV
    tau_YPPV(ii,:) =  (YP(ii)+PV(ii).*SR_range)';
    plot_h(ii,1)=plot(SR_range,tau_YPPV(ii,:),':');
    % PL
    tau_PL_P(ii,:) = (K_P(ii)*SR_range.^(n_P(ii)))';
    set(gca,'ColorOrderIndex',COI);
    plot_h(ii,2)=plot(SR_range,tau_PL_P(ii,:),'--');
    % Reference values at 300 and 600 Fann rpm
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(2),tau_Fann(ii,2),'o');
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(1),tau_Fann(ii,1),'o');
    legendentries{ii,1} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),')'];
    legendentries{ii,2} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),...
        ' \rightarrow n = ', num2str(n_P(ii),2), ', K = ', num2str(K_P(ii),3),')'];
    
    fluidtable{ii,4} = num2str(n_P(ii),2);
    fluidtable{ii,5} = num2str(K_P(ii),2);
    
end
% Plot PAC
if strcmp(PAC,'Y')
    for ii=1:length(fluidlist)
        % Plot Cross fit
        eta = Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*SR_range).^Cross.n{ii});
        plot(SR_range,eta.*SR_range,'-','Color','k');
        % Plot Carreau fit
%         plot(SR,Carreau.mu_inf{ii}+(Carreau.mu_0{ii}-Carreau.mu_inf{ii}).*((1+(Carreau.lambda{ii}.*SR_range).^2).^((Carreau.n{ii}-1)./2)),'-','Color','k');
        pos = max(xlim);
        text(pos-10,...
            Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*(pos+100)).^Cross.n{ii})*(pos+100),fluidlist{ii},...
            'HorizontalAlignment','Right',...
            'VerticalAlignment','Bottom',...
            'Rotation',5);
    end
end
plot_h = [plot_h(:,1); plot_h(:,2)];
legendentries = [legendentries(:,1); legendentries(:,2)];
% legend_h = legend(plot_h,legendentries,'Location','south','FontSize',10);
% % legend_h = legend(plot_h(:,2),legendentries(:,3),'Location','southeast','FontSize',10);
% text(420,max(ylim)-10,'YP [lbf/100ft²]\newlinePV [mPa.s]\newlinen [-]\newlineK[mPa.s^n]');
% Plot text
% colors = get(groot,'defaultAxesColorOrder');
tx = 0.5*(SR_Fann(1)+SR_Fann(2));
ty = 80;
text(tx,ty,...
    {'Typical pipe flow','shear rate range'},...
    'Color','w',...
    'FontSize',18,...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top',...
    'Rotation',0);
% Second set of axis
ax1 = gca; % current axes
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top','YAxisLocation','right',...
    'Xlim',xlim/1.703,'Ylim',ylim/0.4788026,...
    'Xtick',[3 6 100 200 300 600],...
    'Color','none');
xlabel('Fann viscometer speed [rpm]'); ylabel('Shear stress [lbf/100ft²]');
% Minimize white space of figure 
TightFigure([ax1 ax2]); % Handles to axis, if only one just provide one


% Get positions
pos_ax1=get(ax1,'position'); 
pos_ax2=get(ax2,'position'); 
pos_fig=get(gcf,'position');
% % pos_leg=get(legend_h,'position');

% Resize factor
factor = 1 % + pos_leg(3);

% Resize axis
pos_ax1(1)=pos_ax1(1)/factor; pos_ax1(3)=pos_ax1(3)/factor-0.015;
pos_ax2(1)=pos_ax2(1)/factor; pos_ax2(3)=pos_ax2(3)/factor-0.015;

% Resize figure
pos_fig(3)=factor*pos_fig(3)+0.015;
set(gcf,'position',pos_fig);

% Reset axis
set(ax1,'position',pos_ax1);
set(ax2,'position',pos_ax2);

% Programatically move the Legend
% % pos_leg(3)=% % pos_leg(3)/factor; % % pos_leg(1)=1-% % pos_leg(3)-0.015;
% set(legend_h,'Position',% % pos_leg);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'

%% Apparent viscosity vs. shear rate, PL fit of lower shear rate range
PH_FC(3) = CreateFigure( 'Apparent viscosity as a function shear rate w/ PL fit based on annular shear rate range',...
    'Shear rate [1/s]', 'Apparent viscosity [Pa.s]', {'log' 'log'},'PowerPoint');
set(gca,'ColorOrder', ColorSet,'xlim',[1e0 1200], 'ylim',[9e-4 1e1]);
% Shaded shear rate intervalls for pipe and annular flow
SRI = [area([SR_Fann(4) SR_Fann(6)],[max(ylim) max(ylim)])...
area([SR_Fann(2) SR_Fann(1)],[max(ylim) max(ylim)])];
    for aa = 1:length(SRI)
        SRI(aa).FaceColor = [0.8 0.8 0.8];
        SRI(aa).EdgeColor = 'none';
        SRI(aa).FaceAlpha = 0.5;
    end
    
legendentries=cell(length(YP),2);
plot_h=zeros(length(YP),2);
% Plot apparent viscosity of representative drilling fluids
for ii=1:length(YP)
    COI = get(gca,'ColorOrderIndex');
    % YP/PV
    plot_h(ii,1)=plot(SR_range,tau_YPPV(ii,:)./SR_range,':');
    % PL
    set(gca,'ColorOrderIndex',COI);
    plot_h(ii,2)=plot(SR_range,tau_PL_A(ii,:)./SR_range,'--');
    % Reference values at 300 and 600 Fann rpm
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(6),tau_Fann(ii,6)./SR_Fann(6),'o');
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(4),tau_Fann(ii,4)./SR_Fann(4),'o');
    legendentries{ii,1} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),')'];
    legendentries{ii,2} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),...
        ' \rightarrow n = ', num2str(n_A(ii),2), ', K = ', num2str(K_A(ii),3),')'];
end
% Plot PAC
if strcmp(PAC,'Y')
    for ii=1:length(fluidlist)
        % Plot Cross fit
        eta = Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*SR_range).^Cross.n{ii});
        plot(SR_range,eta,'-','Color','k');
        % Plot Carreau fit
%         plot(SR,Carreau.mu_inf{ii}+(Carreau.mu_0{ii}-Carreau.mu_inf{ii}).*((1+(Carreau.lambda{ii}.*SR_range).^2).^((Carreau.n{ii}-1)./2)),'-','Color','k');
        pos = min(xlim)+10^(log10(min(xlim))-1);
        text(pos,....
            Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*pos).^Cross.n{ii}),fluidlist{ii},...
            'HorizontalAlignment','Left',...
            'VerticalAlignment','Bottom',...
            'Rotation',-3);
    end
end
plot_h = [plot_h(:,1); plot_h(:,2)];
legendentries = [legendentries(:,1); legendentries(:,2)];
% legend_h = legend(plot_h,legendentries,'Location','south','FontSize',10);
% % legend_h = legend(plot_h(:,2),legendentries(:,3),'Location','southeast','FontSize',10);
% Plot text
tx = 30;
ty = max(ylim);
txt_A = text(tx,ty,...
    {'Typical annular flow','shear rate range'},...
    'Color','w',...
    'FontSize',18,...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top',...
    'Rotation',0);
% Second set of axis
ax1 = gca; % current axes
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top','YAxisLocation','right',...
    'Xlim',xlim/1.703,'Ylim',ylim*1000,...
    'XScale','log','YScale','log',...
    'Xtick',[3 6 100 200 300 600],...
    'Color','none');
xlabel('Fann viscometer speed [rpm]'); ylabel('Apparent viscosity [cP]');
% Minimize white space of figure 
TightFigure([ax1 ax2]); % Handles to axis, if only one just provide one

% Get positions
pos_ax1=get(ax1,'position'); 
pos_ax2=get(ax2,'position'); 
pos_fig=get(gcf,'position');
% % pos_leg=get(legend_h,'position');

% Resize factor
factor = 1% + pos_leg(3);

% Resize axis
pos_ax1(1)=pos_ax1(1)/factor; pos_ax1(3)=pos_ax1(3)/factor-0.015;
pos_ax2(1)=pos_ax2(1)/factor; pos_ax2(3)=pos_ax2(3)/factor-0.015;

% Resize figure
pos_fig(3)=factor*pos_fig(3)+0.015;
set(gcf,'position',pos_fig);

% Reset axis
set(ax1,'position',pos_ax1);
set(ax2,'position',pos_ax2);

% Programatically move the Legend
% % pos_leg(3)=% % pos_leg(3)/factor; % % pos_leg(1)=1-% % pos_leg(3)-0.015;
% set(legend_h,'Position',% % pos_leg);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'

%% Apparent viscosity vs. shear rate, PL fit of higher shear rate range
PH_FC(4) = CreateFigure( 'Apparent viscosity as a function shear rate w/ PL fit based on pipe shear rate range',...
    'Shear rate [1/s]', 'Apparent viscosity [Pa.s]', {'lin' 'lin'},'PowerPoint');
set(gca,'ColorOrder', ColorSet,'xlim',[0 1200], 'ylim',[0 1e-1]);
% Shaded shear rate intervalls for pipe and annular flow
SRI = [area([SR_Fann(4) SR_Fann(6)],[max(ylim) max(ylim)])...
area([SR_Fann(2) SR_Fann(1)],[max(ylim) max(ylim)])];
    for aa = 1:length(SRI)
        SRI(aa).FaceColor = [0.8 0.8 0.8];
        SRI(aa).EdgeColor = 'none';
        SRI(aa).FaceAlpha = 0.5;
    end
    
legendentries=cell(length(YP),2);
plot_h=zeros(length(YP),2);
% Plot apparent viscosity of representative drilling fluids
for ii=1:length(YP)
    COI = get(gca,'ColorOrderIndex');
    % YP/PV
    plot_h(ii,1)=plot(SR_range,tau_YPPV(ii,:)./SR_range,':');
    % PL
    set(gca,'ColorOrderIndex',COI);
    plot_h(ii,2)=plot(SR_range,tau_PL_P(ii,:)./SR_range,'--');
    % Reference values at 300 and 600 Fann rpm
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(2),tau_Fann(ii,2)./SR_Fann(2),'o');
    set(gca,'ColorOrderIndex',COI);
    plot(SR_Fann(1),tau_Fann(ii,1   )./SR_Fann(1),'o');
    legendentries{ii,1} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),')'];
    legendentries{ii,2} = ['Fluid ' num2str(ii), ' (YP = ', num2str(YP(ii),2), ' , PV = ', num2str(PV(ii)),...
        ' \rightarrow n = ', num2str(n_P(ii),2), ', K = ', num2str(K_P(ii),3),')'];
end
% Plot PAC
if strcmp(PAC,'Y')
    for ii=1:length(fluidlist)
        % Plot Cross fit
        eta = Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*SR_range).^Cross.n{ii});
        plot(SR_range,eta,'-','Color','k');
        % Plot Carreau fit
%         plot(SR,Carreau.mu_inf{ii}+(Carreau.mu_0{ii}-Carreau.mu_inf{ii}).*((1+(Carreau.lambda{ii}.*SR_range).^2).^((Carreau.n{ii}-1)./2)),'-','Color','k');
        pos = max(xlim)-10;
        text(pos,...
            Cross.mu_inf{ii}+(Cross.mu_0{ii}-Cross.mu_inf{ii})./(1+(Cross.lambda{ii}.*pos).^Cross.n{ii}),fluidlist{ii},...
            'HorizontalAlignment','Right',...
            'VerticalAlignment','Bottom',...
            'Rotation',-5);
    end
end
plot_h = [plot_h(:,1); plot_h(:,2)];
legendentries = [legendentries(:,1); legendentries(:,2)];
% legend_h = legend(plot_h,legendentries,'Location','south','FontSize',10);
% % legend_h = legend(plot_h(:,2),legendentries(:,3),'Location','southeast','FontSize',10);
% Plot text
tx = SR_Fann(1)-15; % 0.5*(SR_Fann(1)+SR_Fann(2));
ty = 0.099;
txt_P = text(tx,ty,...
    {'Typical pipe flow','shear rate','range'},...
    'Color','w',...
    'FontSize',18,...
    'FontWeight','normal',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'Rotation',0);
% Second set of axis
ax1 = gca; % current axes
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top','YAxisLocation','right',...
    'Xlim',xlim/1.703,'Ylim',ylim*1000,...
    'Xtick',[3 6 100 200 300 600],...
    'Color','none');
xlabel('Fann viscometer speed [rpm]'); ylabel('Apparent viscosity [cP]');
% Minimize white space of figure 
TightFigure([ax1 ax2]); % Handles to axis, if only one just provide one

% Get positions
pos_ax1=get(ax1,'position'); 
pos_ax2=get(ax2,'position'); 
pos_fig=get(gcf,'position');
% % pos_leg=get(legend_h,'position');

% Resize factor
factor = 1% + pos_leg(3);

% Resize axis
pos_ax1(1)=pos_ax1(1)/factor; pos_ax1(3)=pos_ax1(3)/factor-0.015;
pos_ax2(1)=pos_ax2(1)/factor; pos_ax2(3)=pos_ax2(3)/factor-0.015;

% Resize figure
pos_fig(3)=factor*pos_fig(3)+0.015;
set(gcf,'position',pos_fig);

% Reset axis
set(ax1,'position',pos_ax1);
set(ax2,'position',pos_ax2);

% Programatically move the Legend
% % pos_leg(3)=% % pos_leg(3)/factor; % % pos_leg(1)=1-% % pos_leg(3)-0.015;
% set(legend_h,'Position',% % pos_leg);


% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Produce plot to get legend colors
% 
% CreateFigure( 'Legend colors',...
%     'Shear rate [1/s]', 'Apparent viscosity [Pa.s]', {'lin' 'lin'},'PowerPoint');
% set(gca,'ColorOrder', ColorSet,'xlim',[0 2], 'ylim',[0 11]);
% box off
% grid off
% legendentries=cell(length(YP),2);
% 
% for ii=1:length(YP)
%     COI = get(gca,'ColorOrderIndex');
%     % YP/PV
%     plot(1,11-ii,'o','MarkerFaceColor',ColorSet(ii,:),'MarkerSize',15);
% end

