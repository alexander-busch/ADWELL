%% Get data

A_ScatterAnalysis
% B_Literature


%% Create figure

close all;
fig = figure; % ('units','normalized','outerposition',[0 0 1 1]);
hold on;

%% Formating

xlabel(cat(2,Headers{1,3},' ',Headers{2,3}));
ylabel(cat(2,Headers{1,4},' ',Headers{2,4}));
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [1e-2 1.2e3],...
    'ylim', [2e-2 3e-1],...
    'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1],...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');


% Plot all flow curves

hold on;
for i = 1:13
    plot(PAC4flowcurves(:,3,i),PAC4flowcurves(:,4,i),'ko');
end

% Scatter of Milads data
load Lab2_Scatter2016;
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)+ 3*Lab2_Scatter2016(:,3),'r--');
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2)- 3*Lab2_Scatter2016(:,3),'r--');
plot (Lab2_Scatter2016(:,1),Lab2_Scatter2016(:,2),'ro');

% Laboratory 2 - UiS, Stavanger
% load Lab2;

% Plot raw data
% n = 5;
% Lab2_All = [];
% for i=1:n
%    B = Lab2(2*i-1:2*i,1:length(Lab2));
%    Lab2_All = cat(2,Lab2_All,B);  
% end
% plot(Lab2_All(1,:),Lab2_All(2,:),'k+');

% Plot mean
% MS = 8; % MarkerSize
% MFC = 'k'; % MarkerFaceColor
% B=Lab2(2:2:end,:);
% B=B(4:5,:);
% B_Mean = mean(B,1);
% plot(Lab2(7,:),B_Mean,...
%     'ko','markersize',MS,'MarkerFaceColor',MFC);






%% Create figure

close all;
fig = figure; % ('units','normalized','outerposition',[0 0 1 1]);
hold on;



%% Formating

xlabel(cat(2,Headers{1,3},' ',Headers{2,3}));
ylabel(cat(2,Headers{1,4},' ',Headers{2,4}));
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


%% Data


% plot(Up_SR,Up_AV_M_plus3D,'k--');
plot(Up_SR,Up_AV_M,'ko');
% plot(Up_SR,Up_AV_M_minus3D,'k--');

% plot(Down_SR,Down_AV_M_plus3D,'k--');
plot(Down_SR,Down_AV_M,'ko');
% plot(Down_SR,Down_AV_M_minus3D,'k--');


%% Compute total mean

A_mean = [Up_AV_M, fliplr(Down_AV_M')'];
AV_M = mean(A_mean,2);
plot(Up_SR,AV_M,'ko');


%% Fits

% % createCROSSFits.m;
% 
% % Total mean
% 
%        k =     0.04119
%        n =      0.4662
%        ny_0 =      0.2058
%        ny_inf =   1.13e-08
% 
% plot(Up_SR,ny_inf+(ny_0-ny_inf)./(1+(k.*Up_SR).^n),'k-');
% 
% % Upward sweep
% 
%        k =     0.02878  
%        n =      0.6097  
%        ny_0 =      0.1982  
%        ny_inf =    0.005819  
% 
% plot(Up_SR,ny_inf+(ny_0-ny_inf)./(1+(k.*Up_SR).^n),'k-');
% 
% % Downwards sweep
% 
%        k =      0.1064  
%        n =      0.3598  
%        ny_0 =      0.2254  
%        ny_inf =   3.404e-10  
% 
% plot(Up_SR,ny_inf+(ny_0-ny_inf)./(1+(k.*Up_SR).^n),'k-');

%% Carreau-Yasuda

% Total mean

       a =      0.3541
       k =    0.001275
       n =    0.002558
       ny_0 =       0.214
       ny_inf =       0.001
       
plot(Up_SR,ny_inf+(ny_0-ny_inf)./((1+(k.*Up_SR).^a).^((1-n)./a)),'k-','LineWidth',2);

% Upward sweep

       a =      0.5785 
       k =     0.02459  
       n =      0.4042  
       ny_0 =      0.1992 
       ny_inf =       0.001 

plot(Up_SR,ny_inf+(ny_0-ny_inf)./((1+(k.*Up_SR).^a).^((1-n)./a)),'k-');

% Downwards sweep

       a =      0.2462 
       k =   0.0001924 
       n =    0.006498  
       ny_0 =      0.2455  
       ny_inf =       0.001  

plot(Up_SR,ny_inf+(ny_0-ny_inf)./((1+(k.*Up_SR).^a).^((1-n)./a)),'k-');





%% 













































