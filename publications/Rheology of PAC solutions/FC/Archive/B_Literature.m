

% n=8;
% figure();
% hold on;
% for i=1:n
%     plot(ShearRate,PowerLawCoefficients{i,2}.*ShearRate.^(PowerLawCoefficients{i,3}-1));
% end
% legend('Upward Mean', 'Upward Mean + 3SD', 'Upward Mean - 3SD', 'Downward Mean', 'Downward Mean + 3SD', 'Downward Mean - 3SD', PowerLawCoefficients{:,1},'Location','southwest');
% set(gca,'XScale','log','YScale','log');



%% Create figure

close all;
fig = figure; % ('units','normalized','outerposition',[0 0 1 1]);
hold on;


%% Plot this study results 

% Laboratory 1 - SINTEF Petroleum AS, Bergen
Lab1 = zeros(2,2*length(Up_SR));
Lab1(1,:) = cat(2,Up_SR.',Down_SR.');
Lab1(2,:) = cat(2,Up_AV_M_plus3D.',Down_AV_M_minus3D.');
plot(Lab1(1,:),Lab1(2,:),'k--');

% load Fann;
% plot (Fann(1,:),Fann(2,:)./Fann(1,:)+ 3*Fann(3,:)./Fann(1,:),'r--');
% plot (Fann(1,:),Fann(2,:)./Fann(1,:)- 3*Fann(3,:)./Fann(1,:),'r--');



% Laboratory 2 - UiS, Stavanger
load Lab2_Scatter.mat;
% Plot raw data
% n = 5;
% Lab2_All = [];
% for i=1:n
%    B = Lab2(2*i-1:2*i,1:length(Lab2));
%    Lab2_All = cat(2,Lab2_All,B);  
% end
% plot(Lab2_All(1,:),Lab2_All(2,:),'k+');
% Plot mean
MS = 8; % MarkerSize
MFC = 'k'; % MarkerFaceColor
B=Lab2_Scatter(2:2:end,:);
B=B(4:5,:);
B_Mean = mean(B,1);
plot(Lab2_Scatter(7,:),B_Mean,...
    'ko','markersize',MS,'MarkerFaceColor',MFC);

% load Lab2_Scatter2016; % This seems to be rather Time (2015) data...
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)+ 3*Lab2_Scatter2016(:,3),'r--');
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2),'r--');
% plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)- 3*Lab2_Scatter2016(:,3),'r--');

% load Lab3;
% plot (Lab3(:,1),Lab3(:,2));


%% Plot literature benchmark

addpath C:\Users\alexabus\IEPT1325_Local\Data\SW_working_directories\MATLAB\AdWell\Rheology\PAC\FC\Literature_PAC4;


LW = 1; % LineWidth
MS = 8; % MarkerSize
MFC = 'w'; % MarkerFaceColor

% Time (2015)
load Time2015_CC;
plot (Time2015_CC(:,1),Time2015_CC(:,2),...
    'ko','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);
% load Time2015_CP;
% plot (Time2015_CP(:,1),Time2015_CP(:,2),...
%     'r+','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Johnsen (2014)
load Johnsen2014;
plot(Johnsen2014(1,:),Johnsen2014(2,:),...
    'k^','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Time and Rabenjafimanantsoa (2012)
load Time2012;
plot(Time2012(1,:),Time2012(2,:),...
    'ks','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);

% Rabenjamanantsoa et al. (2011)
load Rabenjamanantsoa2011;
plot(Rabenjamanantsoa2011(1,:),Rabenjamanantsoa2011(2,:),...
    'kd','MarkerSize',MS,'MarkerFaceColor',MFC,'LineWidth',LW);


%% Formating

xlabel(cat(2,headers{1,3},' ',headers{2,3}));
ylabel(cat(2,headers{1,4},' ',headers{2,4}));
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [0 1200],...
    'ylim', [2e-2 3e-1],...
    'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1],...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');


%% Legend

h_legend = legend(...
    'This study: Lab1 (Up- & Downward Mean \pm 3SD)', ...
    'This study: Lab2 (Mean)', ...
    'Time et al. (2015)', ...
    'Johnsen (2014)', ...
    'Time et al. (2012)', ...
    'Rabenjamanantsoa et al. (2011)');

set(h_legend,'Location','southwest','FontSize',22);




%% Print as image


screen_size = get(0, 'ScreenSize');
origSize = get(fig, 'Position'); % grab original on screen size
% set(fig, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
set(fig,'Units','centimeters','Position',[0 0 32 18])
 
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% print(fig,'M:\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper\Figures\Figure2_Flowcurves','-depsc');  % '-dpng' '-depsc' 'jpeg' 'eps'
    
%  saveas(eps, 'M:\data\AdWell\8 - Publications\2016-05 - Rheological properties of PAC\Paper\Figures\Figure2_Flowcurves') % save figure
% set(fig,'Position', origSize) %set back to original dimensions




% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
%
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 26 14])
% 
% print(fig,'M:\data\AdWell\8 - Publications\2016-05 - Rheological properties of PAC\Paper\Figures\Figure2_Flowcurves','-depsc');  % '-dpng' '-depsc' 'jpeg' 'eps'


