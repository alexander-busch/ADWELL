%% Create figure
A_ScatterAnalysis;

close all;
fig = figure('Units','centimeters','Position',[1 1 21 14.8]); % DIN A5 
hold on;

%% Formating

% xlabel(cat(2,headers{1,3},' ',headers{2,3}));
% ylabel(cat(2,headers{1,4},' ',headers{2,4}));
xlabel('$\dot\gamma$ [1/s]','FontSize',12,'Interpreter','latex');
ylabel('$\eta [Pa.s]$','FontSize',12,'Interpreter','latex');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [0.01 1200],...
    'ylim', [9e-3 2e-0],... %3e-1
    'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    'box','on',...
    'FontSize',12);
    % 'FontSize',24);
set(gcf,...
    'color','w');

col_PAC8 = [0,0,0]; %[0,0,1]; 
col_PAC4 = [0,80,158]/255; %[0,0.5,1]; 
col_PAC2 = [144, 73, 45]/255; %[0,1,1]; 


%% PAC2 this study

% PAC2 results are taken from UiS as Bergen data is wrong (PAC8 instead of
% PAC2)

% PAC 2 - Laboratory 2 - UiS - insufficient shear rate range
load pac2_1;
load pac2_2;
pac2(1,:) = pac2_1(:,2);
pac2(2,:) = pac2_2(:,2);
%pac2 = mean(pac2,2);
pac2 = (pac2(1,:)+pac2(2,:))./2;
Lab2_PAC2_M_Up = plot(pac2_1(1:25),pac2(1:25)','-','Color',col_PAC2);
Lab2_PAC2_M_Down = plot(pac2_1(26:50),pac2(26:50)','-','Color',col_PAC2);
ShearRate_PAC2=pac2_1(1:25);
ApparentViscosity_PAC2=(pac2(1:25)'+flip(pac2(26:50)'))./2;

%% PAC2 literature benchmark

LW = 1; % LineWidth
MS = 6; % MarkerSize
MFC = 'w'; % MarkerFaceColor

% Khatibi and Time (2016)
load Khatibi2016_PAC2;
Khatibi2016_PAC2 = plot (Khatibi2016_PAC2(:,1),Khatibi2016_PAC2(:,2)./1000,...
    'o','Color',col_PAC2,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Zoric et al. (2015) - Equivalent to Time et al. (2012)
% load Zoric2015_PAC2;
% plot (Zoric2015_PAC2(:,1),Zoric2015_PAC2(:,2),...
%      'go','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time and Rabenjafimanantsoa (2012) - A - Equivalent to Time et al. (2012)
load Time2012A_PAC2;
Time2012A_PAC2 = plot(Time2012A_PAC2(:,1),Time2012A_PAC2(:,2),...
    '^','Color',col_PAC2,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time and Rabenjafimanantsoa (2012) - S
load Time2012S_PAC2;
Time2012S_PAC2 = plot(Time2012S_PAC2(:,1),Time2012S_PAC2(:,2),...
    's','Color',col_PAC2,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Rabenjamanantsoa et al. (2011) - Do not plot as this leads to size issues
% with the legend
% load Rabenjamanantsoa2011_PAC2;
% Rabenjamanantsoa2011_PAC2 = plot(Rabenjamanantsoa2011_PAC2(:,1),Rabenjamanantsoa2011_PAC2(:,2),...
%     'kd','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);




% Sigve 2018 NTNU Petroleum lab

%PAC2 = importPAC2data('E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\PAC\Resub1\FC\Data\Lab3_NTNU\180205_PAC2_FC_Sigve.xlsx','Ark1','AJ8:AL29');
%PAC2_M = plot(PAC2.SR,PAC2.AV,'r--','LineWidth',2);


%% PAC4 & PAC8 this study results


% PAC 4 - Laboratory 1 - SINTEF Petroleum AS, Bergen
% Lab1 = zeros(2,2*length(Up_SR));
% Lab1(1,:) = cat(2,Up_SR.',Down_SR.');
% Lab1(2,:) = cat(2,Up_AV_M.',Down_AV_M.');
% plot(Lab1(1,:),Lab1(2,:),'k--');

Lab1_PAC4_M_Up = plot(Up_SR,PAC4_AV_M(1:50),'-','Color',col_PAC4); % Upward mean
Lab1_PAC4_M_Down = plot(Down_SR,PAC4_AV_M(61:110),'-','Color',col_PAC4); % Downward mean
Lab1_PAC4_M_Up3SD = plot(Up_SR,PAC4_AV_M(1:50)+3.*PAC4_AV_SD(1:50),':','Color',col_PAC4); % Upward mean + 3SD
Lab1_PAC4_M_Down3SD = plot(Down_SR,PAC4_AV_M(61:110)-3.*PAC4_AV_SD(61:110),':','Color',col_PAC4); % Downward mean - 3SD


% Viscosity ratio: Downward / Upward
AppVisRatio = PAC4_AV_M(61:110)'./fliplr(PAC4_AV_M(1:50)');


% Uncertainty definition, usuallz standard uncertaintz u = 1SD
% Here we use 3SD
factor = 3;

% Relative SD of upward sweep and downward sweeps and mean relative SD
Lab1_SDr_up = factor*PAC4_AV_SD(1:50)./PAC4_AV_M(1:50);
Lab1_SDr_down = factor*PAC4_AV_SD(61:110)./PAC4_AV_M(61:110);
Lab1_SDr = [Lab1_SDr_up, fliplr(Lab1_SDr_down')'];
Lab1_SDr = mean(Lab1_SDr ,2);

% Standard error of upward and downward sweeps
Lab1_SE_up = PAC4_AV_SD(1:50)./sqrt(13);
Lab1_SE_down = PAC4_AV_SD(61:110)./sqrt(13);

% Relative standard error of upward and downward sweeps
Lab1_SEr_up = Lab1_SE_up./PAC4_AV_M(1:50);
Lab1_SEr_down = Lab1_SE_down./PAC4_AV_M(61:110);

% Total mean of upward and downward sweep
Lab1_PAC4_M = [PAC4_AV_M(1:50), fliplr(PAC4_AV_M(61:110)')'];
Lab1_PAC4_M = mean(Lab1_PAC4_M,2);

% "Standard deviation" for total mean of upward and downward sweep
Lab1_SD_up_M = factor*PAC4_AV_SD(1:50)+(PAC4_AV_M(1:50)-Lab1_PAC4_M);
Lab1_SD_down_M = fliplr(factor*PAC4_AV_SD(61:110)')'+(Lab1_PAC4_M-fliplr(PAC4_AV_M(61:110)')');

% Realative "Standard deviation" for total mean of upward and downward
% sweep and mean relative SD
Lab1_SDr_up_M = Lab1_SD_up_M./Lab1_PAC4_M;
Lab1_SDr_down_M = Lab1_SD_down_M./Lab1_PAC4_M;
Lab1_SDr_M = [Lab1_SDr_up_M, Lab1_SDr_down_M];
Lab1_SDr_M = mean(Lab1_SDr_M,2);

% Checks
% k = find(Up_SR==140) 
% k = find(Up_SR==33.4) 
% k = find(Up_SR==4.96) 
% PAC4_AV_M(k)
% Lab1_PAC4_M(k)
% check = Lab1_SDr_up_M./Lab1_SDr_up
% check = Lab1_SDr_M./Lab1_SDr 
% check = Lab1_SDr_M
% check(k)

Lab1_PAC4_M = plot(Up_SR,Lab1_PAC4_M,'--','Color',col_PAC4);

% load Fann;
% plot (Fann(1,:),Fann(2,:)./Fann(1,:)+ 3*Fann(3,:)./Fann(1,:),'r--');
% plot (Fann(1,:),Fann(2,:)./Fann(1,:)- 3*Fann(3,:)./Fann(1,:),'r--');


% PAC 8 - Laboratory 1 - SINTEF Petroleum AS, Bergen
Lab1_PAC8_M_Up = plot(Up_SR,PAC8_AV_M(1:50),'-','Color',col_PAC8);
Lab1_PAC8_M_Down = plot(Down_SR,PAC8_AV_M(61:110),'-','Color',col_PAC8);
% % Total mean
% pac8 (1,:) = PAC8_AV_M(1:50)
% pac8 (2,:) = flip(PAC8_AV_M(61:110))
% pac8 = mean(pac8(:,:),1)
% plot(Up_SR,pac8','k--');

% PAC 4 - Laboratory 2 - UiS, Stavanger
load Lab2;
% First five from 150430_Rheo_PAC200&PAC400_Analysed.xlsx
% Identical to 150605_Rheo_PAC200&PAC400_.xlsx



% Plot raw data
% figure();
% hold on;
% grid('on');
% Lab2_All = [];
% for i=1:5
%     B = Lab2(2*i-1:2*i,1:length(Lab2));
%     Lab2_All = cat(2,Lab2_All,B);
% end
% plot(Lab2_All(1,:),Lab2_All(2,:),'k+');
% set(gca,...
%     'XScale','log',...
%     'YScale','log',...
%     'xlim', [0.01 1200],...
%     'ylim', [9e-3 2e-0],... %3e-1
%     'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
%     'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
%     'box','on',...
%     'FontSize',24);

% Plot mean
MFC = 'k'; % MarkerFaceColor
B=Lab2(2:2:end,:);
B=B(4:5,:);
ShearRate = Lab2(7,:);
ApparentViscosity = mean(B,1);
%ApparentViscosity = (B(1,:)+B(2,:)./2);
% Lab2_PAC4_M = plot(ShearRate,ApparentViscosity,'ko','markersize',MS,'MarkerFaceColor',MFC);

ShearRate_PAC4 = ShearRate(1:20);
ApparentViscosity = [ApparentViscosity(1:20); fliplr(ApparentViscosity(21:40))];
ApparentViscosity_PAC4 = mean(ApparentViscosity,1);
Lab2_PAC4_M = plot(ShearRate_PAC4,ApparentViscosity_PAC4,'--','Color',col_PAC4,'LineWidth',2);

% load Lab2_Scatter2016; % This seems to be rather Time (2015) data...
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)+ 3*Lab2_Scatter2016(:,3),'r--');
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2),'r--');
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)- 3*Lab2_Scatter2016(:,3),'r--');
% 

% % PAC 4 - Laboratory 3 - NTNU, Trondheim
% load Lab3;
% plot (Lab3(:,1),Lab3(:,2));

parentpath = cd(cd('..'));
fig_path = [parentpath '\NormalStressDifferences'];
save([fig_path '\FC.mat'],...
    'ShearRate_PAC2',...
    'ApparentViscosity_PAC2',...
    'ShearRate_PAC4',...
    'ApparentViscosity_PAC4');

clear('ShearRate',...
    'ApparentViscosity');


%% PAC4 literature benchmark

MFC = 'w'; % MarkerFaceColor

% Khatibi and Time (2016) - Equivalent to lab2 data
% load Khatibi2016_PAC4;
% plot (Khatibi2016_PAC4(:,1),Khatibi2016_PAC4(:,2)./1000,...
%     'ro','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Zoric et al. (2015) - Equivalent to Time et al. (2012)
% load Zoric2015_PAC4;
% figure
% plot (Zoric2015_PAC4(:,1),Zoric2015_PAC4(:,2),...
%     'go','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time (2015)
load Time2015_PAC4_CC;
Time2015_PAC4_CC = plot (Time2015_PAC4_CC(:,1),Time2015_PAC4_CC(:,2),...
    'o','Color',col_PAC4,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time (2014) - Equivalent to Johnsen (2014)
% load Time2014_PAC4;
% hold on
% plot (Time2014_PAC4(:,1),Time2014_PAC4(:,2)./1000,...
%     'bo','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Johnsen (2014)
load Johnsen2014_PAC4;
Johnsen2014_PAC4 = plot(Johnsen2014_PAC4(:,1),Johnsen2014_PAC4(:,2),...
    '^','Color',col_PAC4,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time and Rabenjafimanantsoa (2012) - A - Equivalent to Time et al. (2012)
% load Time2012A_PAC4;
% plot(Time2012A_PAC4(:,1),Time2012A_PAC4(:,2),...
%     'ks','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time and Rabenjafimanantsoa (2012) - S
load Time2012S_PAC4;
Time2012S_PAC4 = plot(Time2012S_PAC4(:,1),Time2012S_PAC4(:,2),...
    's','Color',col_PAC4,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Rabenjamanantsoa et al. (2011) - Do not plot as this leads to size issues
% with the legend
% load Rabenjamanantsoa2011_PAC4;
% Rabenjamanantsoa2011_PAC4 = plot(Rabenjamanantsoa2011_PAC4(:,1),Rabenjamanantsoa2011_PAC4(:,2),...
%     'kd','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

%% PAC8 literature benchmark

% Time (2014)
load Time2014_PAC8;
hold on
Time2014_PAC8 = plot (Time2014_PAC8(:,1),Time2014_PAC8(:,2)./1000,...
    'o','Color',col_PAC8,'MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);


%% Legend

% Toggle legend display status
set(get(get(Lab2_PAC2_M_Up,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab2_PAC2_M_Down,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC4_M_Up,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC4_M_Down,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC4_M_Up3SD,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC4_M_Down3SD,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab2_PAC4_M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC8_M_Up,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC8_M_Down,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(Lab1_PAC4_M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Display legend

h_legend = legend(...
    'Khatibi et al. (2016)',...
    'Time et al. (2012)', ...
    'Time et al. (2012)', ...
    'Time et al. (2015)', ...
    'Johnsen (2014)', ...
    'Time et al. (2012)', ...
    'Time et al. (2014)');
  %  'PAC2 (Rabenjamanantsoa et al. (2011)',...
  %  'PAC4 (Rabenjamanantsoa et al. (2011)',...
  
set(h_legend,'Location','northeast','FontSize',12); % 'northeastoutside'
legend boxoff;


%% Labeling of flow curves

x1 = 1.1e-2
y1 = 1.25;
txt1 = 'PAC8';
text(x1,y1,txt1,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','bold','Color',col_PAC8); %'interpreter','latex'

y1 = 0.2;
txt1 = 'PAC4';
text(x1,y1,txt1,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'FontWeight','bold','Color',col_PAC4); %'interpreter','latex'

x1 = 0.12;
y1 = 0.051;
txt1 = 'PAC2';
text(x1,y1,txt1,'HorizontalAlignment','right','FontSize',16,'FontWeight','bold','Color',col_PAC2); %'interpreter','latex'


%% Labeling of upward & downward sweeps curves

% PAC8
x1 = 1;
y1 = 1.28;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-10,'Color',col_PAC8); %'interpreter','latex'
annotation('textarrow',[0.36 0.5],[0.9 0.86],'String','','Color',col_PAC8);
y1 = 0.65;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-10,'Color',col_PAC8); %'interpreter','latex'
annotation('textarrow',[0.5 0.36],[0.79 0.83],'String','','Color',col_PAC8);


% PAC4
x1 = 22;
y1 = 0.16;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-21,'Color',col_PAC4); %'interpreter','latex'
annotation('textarrow',[0.6 0.74],[0.61 0.53],'String','','Color',col_PAC4);
y1 = 0.075;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-15,'Color',col_PAC4); %'interpreter','latex'
annotation('textarrow',[0.74 0.6],[0.435 0.505],'String','','Color',col_PAC4);


% PAC2
x1 = 22;
y1 = 0.06;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-15,'Color',col_PAC2); %'interpreter','latex'
annotation('textarrow',[0.6 0.74],[0.44 0.38],'String','','Color',col_PAC2);
y1 = 0.028;
text(x1,y1,'','HorizontalAlignment','center','FontSize',14,'Rotation',-12,'Color',col_PAC2); %'interpreter','latex'
annotation('textarrow',[0.74 0.6],[0.28 0.33],'String','','Color',col_PAC2);





%% Rheometer accuracy
accuracy;











%% Box for legend
% 
% % PAC2 box
% xmin = 1.05e2;
% xmax = 2e2;
% ymax = 1.85;
% ymin = 1.5;
% 
% plot([xmin xmax], [ymax ymax],'k-')
% plot([xmin xmin], [ymin ymax],'k-')
% plot([xmin xmax], [ymin ymin],'k-')
% 
% x1 = 1.14e2;
% y1 = (ymax-ymin)/(log(ymax)-log(ymin));
% txt1 = 'PAC8';
% txtlabel = text(x1,y1,txt1)
% set(txtlabel,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',22,...
%     'FontWeight','bold',...
%     'rotation', 90); %'interpreter','latex'
% 
% 
% % PAC4 box
% xmin = 1.05e2;
% xmax = 2e2;
% ymax = ymin;
% ymin = 2.7e-1;
% 
% plot([xmin xmax], [ymax ymax],'k-')
% plot([xmin xmin], [ymin ymax],'k-')
% plot([xmin xmax], [ymin ymin],'k-')
% 
% x1 = 1.14e2;
% y1 = (ymax-ymin)/(log(ymax)-log(ymin));
% txt1 = 'PAC4';
% txtlabel = text(x1,y1,txt1)
% set(txtlabel,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',22,...
%     'FontWeight','bold',...
%     'rotation', 90); %'interpreter','latex'
% 
% 
% % PAC8 box
% xmin = 1.05e2;
% xmax = 2e2;
% ymax = ymin;
% ymin = 2.19e-1;
% 
% plot([xmin xmax], [ymax ymax],'k-')
% plot([xmin xmin], [ymin ymax],'k-')
% plot([xmin xmax], [ymin ymin],'k-')
% 
% x1 = 1.14e2;
% y1 = (ymax-ymin)/(log(ymax)-log(ymin));
% txt1 = 'PAC8';
% txtlabel = text(x1,y1,txt1)
% set(txtlabel,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',22,...
%     'FontWeight','bold',...
%     'rotation', 90); %'interpreter','latex'

%% Print as image

% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.003;
factor_ver = 0.01;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];

% Specify Figure Size and Page Size
set(fig,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
fig_path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub2\Figures\';
fig_name = 'FC';
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig,[fig_path num2str(3)],'-dpng','-r600') 
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file
