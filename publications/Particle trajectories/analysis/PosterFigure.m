% Clean up
clear all;
close all;
clc;

addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\export_fig');
% addpath('O:\OneDrive_NTNU\SWwd\MATLAB\Generic\export_fig');

% todo = 1; % read EXP data and setup CFD cases
todo = 0; % read CFD data and plot all results

v_ini = 0; % Compute vx0 and vy0 from experimental data
% v_ini = 1; % As given in paper

% Path of Fluent exports folder
% 'E:\... = Laptop
% 'O:\... = Desktop
resultspath = 'E:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - 2D\archive\170914_2D_C_expVelProf_slip_Cross_SR=1\exports';

% Path for figure export
figurepath = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2017-11 - A DPM drag model for shear-thinning fluids';

% Definitions & parameters
Definitions;
Parameters;

% Tables
CFD_asis = table;
CFD_new = table;


% Current row
row = 1;

% Build tables, loop all fluids, particle diameters and bulk velocities
for i = 1:length(fluidlist)
    for j = 1:length(d_p)
        for k = 1:length(U)

            % Build CFD tables
            CFDdata2Table;
  
            % Increase table row index
            row = row + 1;
        end % of k
    end % of j
end % of i


% Delete empty rows = non-existing experiments 
toDelete = CFD_new.delete == 1;
CFD_asis(toDelete,:) = [];
CFD_new(toDelete,:) = [];

save('CFD.mat',...
    'CFD_asis',...
    'CFD_new');



%% Trajectories


load('Khatibi_et_al_2016.mat');
load('CFD.mat');
addpath E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files;
i_max = size(Cases);

fig_titles = {'CFD vs. Exp. results','Case setup and velocity profiles'};
label_x = {'Fluid x-velocity u_x [m/s]'; 'Channel x-coordinate [m]'};
label_y = 'Channel y-coordinate [m]';
fig_fontsize = 30;

fig = figure('Name',fig_titles{1},'color','w','Units','centimeters','Position',[1 1 105.125 25.0]);
hold on;
grid on;
box on;
xlabel(label_x{2},'interpreter','tex','FontName','Arial');
ylabel(label_y,'interpreter','tex','FontName','Arial');

set(gca,...
    'FontSize',fig_fontsize,...
    'xLim',[-0.02 0.12],...
    'yLim',[0 0.02],...
    'xTick',[0 0.02 0.04 0.06 0.08 0.10 0.12]);

plot([-0.02 0.12], [0 0],'k-','LineWidth',6);
plot([-0.02 0.12], [0.02 0.02],'k-','LineWidth',6);


% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.015;
factor_ver = 0.015;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];


% Loop all fluids, bulk velocities and particle diameters and plot
% results
for i = 1:i_max(1) % i = 2, i = 3, i = 4

    % Subplot
    if Cases.fluid(i,1) == 'H2O'           
        plotcolor = colorlist{1};
    elseif Cases.fluid(i,1) == 'PAC2'         
        plotcolor = colorlist{2};
    elseif Cases.fluid(i,1) == 'PAC4'
        plotcolor = colorlist{3};
    else
    end

    % Define markersize (= particle size)
    if Cases.d_p(i) == 0.00116
        markersize = markersizelist(1);
    elseif Cases.d_p(i) == 0.002
        markersize = markersizelist(2);
    else
        markersize = markersizelist(3);
    end

    % Check the length of vectors to be plotted as sometimes Fluent exports
    % inconsistent data in terms of array elements
    length_X = length(CFD_asis.X{i,1});
    length_Y = length(CFD_asis.Y{i,1});

    if length_X < length_Y
        CFD_asis.X{i,1}(length_X+1) = CFD_asis.X{i,1}(length_X);
        CFD_asis.t{i,1}(length_X+1) = CFD_asis.t{i,1}(length_X);    

    elseif length_X > length_Y
        CFD_asis.Y{i,1}(length_Y+1) = CFD_asis.Y{i,1}(length_Y);
    else
    end

    length_X = length(CFD_new.X{i,1});
    length_Y = length(CFD_new.Y{i,1});

    if length_X < length_Y
        CFD_new.X{i,1}(length_X+1) = CFD_new.X{i,1}(length_X);
        CFD_new.t{i,1}(length_X+1) = CFD_new.t{i,1}(length_X);    

    elseif length_X > length_Y
        CFD_new.Y{i,1}(length_Y+1) = CFD_new.Y{i,1}(length_Y);
    else
    end

    if i == 5
        % Skip PAC4, U = 0.048
    elseif i == 7
        % Skip PAC4, U = 0.048
    else    
        plot(CFD_asis.X{i,1}-CFD_asis.X{i,1}(1),CFD_asis.Y{i,1},'Color',plotcolor,'LineStyle','--','LineWidth',2); % 1 is the row of the table and represents a particular case
        plot(CFD_new.X{i,1}-CFD_new.X{i,1}(1),CFD_new.Y{i,1},'Color',plotcolor,'LineStyle','-','LineWidth',2); % 1 is the row of the table and represents a particular case
        % plot(CFD_Khatibi.X{i,1},CFD_Khatibi.Y{i,1},'Color',plotcolor,'LineStyle',':'); % 1 is the row of the table and represents a particular case
        scatter(EXP.X{i,1},EXP.Y{i,1},markersize,'MarkerEdgeColor',plotcolor,'LineWidth',2); % 1 is the row of the table and represents a particular case
    end

end % of i


% Inlet Velocity profile
H = 0.02; % channel_height
h = 0.0; % bed height
U_l = 0.008; % fluid bulk vel
y = linspace(0,H);

u=2.*U_l.*(1-(y-H./2).^2./(H./2).^2);  % h = channel y-coordinate from 0...H , H = channel height
plot(u-0.02,y,'color','[0.00392 0.31373 0.6196]');

% Arrows
for i = 1:19 % i = 20
    y = i/1000;
    u=2.*U_l.*(1-(y-H./2).^2./(H./2).^2);
    arrow([-0.02 y],[u-0.02  y],'color','[0.00392 0.31373 0.6196]');
end


% Save figures
savefig(fig, [figurepath '\ParticleTrajectories.fig']);
% set(gcf,'PaperPosition',[1 1 105.125 25.0]);
% set(gcf,'PaperSize',[105.125 25.0]);
set(gcf,'PaperPositionMode','auto'); % prints size displayed on screen
print(fig,[figurepath '\ParticleTrajectories'],'-depsc2');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[figurepath '\ParticleTrajectories'],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
set(gcf, 'color', 'none');
set(findall(gcf,'type','axes'), 'color', 'none');
export_fig ParticleTrajectories.png;
movefile('ParticleTrajectories.png', figurepath);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
    
    
%% Plot exp velocity profiles at inlet

fig_fontsize = 20;
rows = [1, 3, 6];

x0 = 0.02;
y0 = 0.00;

fig = figure('Name',fig_titles{1},'color','w','Units','centimeters','Position',[1 1 39 23.6]);
hold on;
grid on;
box on;
xlabel(label_x{1},'interpreter','tex');
ylabel(label_y,'interpreter','tex');
set(gca,...
    'FontSize',fig_fontsize,...
    'xLim',[0 0.16],...
    'yLim',[0 0.02]);

% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.015;
factor_ver = 0.015;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];


fig_fontsize = 30;

for i=1:length(rows)

    x = EXP.y{rows(i),1} - x0;
    y = EXP.u{rows(i),1};
        
    if i == 1
        Fit = CreateVelocityProfileFit_H2O(x,y);
    else
        Fit = CreateVelocityProfileFit_PAC(x,y);
    end
    
    scatter(y,x+x0,'MarkerEdgeColor',colorlist{i},'Marker','o','LineWidth',2);
    x = linspace(-x0,0);
   
    if i == 1
        plot(Fit.p1.*x.^5 + Fit.p2.*x.^4 + Fit.p3.*x.^3 + Fit.p4.*x.^2 + Fit.p5.*x + Fit.p6,x+x0,'Color',colorlist{i},'LineStyle','-','LineWidth',2);
        str = {'H2O, U = 0.085 m/s'};% ,...
            %['u_x(y) = ' num2str(Fit.p1,2) ' Y^5 ' num2str(Fit.p2,2) ' Y^4 ' num2str(Fit.p3,2) ' Y^3 ' num2str(Fit.p4,2) ' Y^2 ' num2str(Fit.p5,2) ' Y +' num2str(Fit.p6,2)]};
        text(0.065,0.0188,str,'Color','blue','HorizontalAlignment' ,'left','VerticalAlignment','Bottom','FontSize',fig_fontsize,'interpreter','tex');
    else
        plot(Fit.p1.*x.^2 + Fit.p2.*x + Fit.p3,x+x0,'Color',colorlist{i},'LineStyle','-','LineWidth',2);
        if i == 2
            str = {'PAC2, U = 0.048 m/s'};% ,...
                %['u_x(y) = ' num2str(Fit.p1,2) ' Y^2 ' num2str(Fit.p2,2) ' Y +' num2str(Fit.p3,2)]};
            text(0.06,0.0038,str,'Color','green','HorizontalAlignment' ,'right','VerticalAlignment','Bottom','FontSize',fig_fontsize,'interpreter','tex');
        else
            str = {'PAC4, U = 0.085 m/s'};% ,...
                % ['u_x(y) = ' num2str(Fit.p1,2) ' Y^2 ' num2str(Fit.p2,2) ' Y +' num2str(Fit.p3,2)]};
            text(0.116,0.0152,str,'Color','red','HorizontalAlignment' ,'left','VerticalAlignment','Bottom','FontSize',fig_fontsize,'interpreter','tex');
        end
        

    end

    
end


% Save figures
savefig(fig, [figurepath '\VelocityProfiles.fig']);,
set(gcf,'PaperPositionMode','auto'); % prints size displayed on screen
print(fig,[figurepath '\VelocityProfiles'],'-depsc2');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[figurepath '\VelocityProfiles'],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
set(gcf, 'color', 'none');
set(findall(gcf,'type','axes'), 'color', 'none');
export_fig VelocityProfiles.png;
movefile('VelocityProfiles.png', figurepath);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');


%% Rheology
fig_fontsize = 20;
rows = [1, 3, 6];

x0 = 0.02;
y0 = 0.00;

fig = figure('Name','Rheometric data','color','w','Units','centimeters','Position',[1 1 39 18]);
hold on;
grid on;
box on;
xlabel('Shear rate \gamma'' [1/s]','interpreter','tex');
ylabel('Apparent viscosity \eta_f [Pa.s]','interpreter','tex');
set(gca,...
    'FontSize',fig_fontsize,...
    'Xscale', 'log',...
    'Yscale', 'log',...
    'Ylim', [0 1e0]);

% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.015;
factor_ver = 0.015;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];

% Create legend cell array
Legend = cell(2*length(fluidlist),1);

% Create Shear rate vector for plotting Cross fit
SR = logspace(-2,4);

% Load table for model coefficients
load Cross;
load Carreau;
load RheoData;


% Loop all fluids, get rheometric data, fit model, plot
for i = 1:length(fluidlist)

    % Color
%     if Cases.fluid(i,1) == 'H2O'           
%         plotcolor = colorlist{1};
%     elseif Cases.fluid(i,1) == 'PAC2'         
%         plotcolor = colorlist{2};
%     elseif Cases.fluid(i,1) == 'PAC4'
%         plotcolor = colorlist{3};
%     else
%     end
    
    % Plot rheometric data
    if i == 2
        a = 3;
    elseif i == 3
        a = 5;
    else
        a = 1;
    end
    scatter(RheoData.gamma_dot{a}, RheoData.app_vis{a},'o','MarkerEdgeColor',colorlist{i},'LineWidth',2);
    
    % Plot Cross fit
    plot(SR,Cross.mu_inf{i}+(Cross.mu_0{i}-Cross.mu_inf{i})./(1+(Cross.lambda{i}.*SR).^Cross.n{i}),'-','Color',colorlist{i},'LineWidth',2);

    % Plot Carreau fit
    % plot(SR,Carreau.mu_inf{i}+(Carreau.mu_0{i}-Carreau.mu_inf{i}).*((1+(Carreau.lambda{i}.*SR).^2).^((Carreau.n{i}-1)./2)),'-','Color',colorlist{i},'LineWidth',2);
  
    % Create legend element
    if i == 1
    else
        
    end
    Legend{2*i-1} = [fluidlist{i} '; Rheo. data of Khatibi et al. (2016)'];
    Legend{2*i-0} = [fluidlist{i} '; Cross fit (\lambda = ',...
        num2str(Cross.lambda{i},3) ', n = ',...
        num2str(Cross.n{i},3) ', \mu_0 = ',...
        num2str(Cross.mu_0{i},3) ', \mu_i_n_f = ',...
        num2str(Cross.mu_inf{i},3) ')'];
%     Legend{3*i} = [fluidlist{i} '; Carreau fit (\lambda = ',...
%         num2str(Carreau.lambda{i},3) ', n = ',...
%         num2str(Carreau.n{i},3) ', \mu_0 = ',...
%         num2str(Carreau.mu_0{i},3) ', \mu_i_n_f = ',...
%         num2str(Carreau.mu_inf{i},3) ')'];
end

% Display legend
Legend = legend(Legend);
set(Legend,...
    'Location','southwest',...
    'FontSize',20,...
    'interpreter','tex');


% Save figures
savefig(fig, [figurepath '\RheometricData.fig']);,
set(gcf,'PaperPositionMode','auto'); % prints size displayed on screen
print(fig,[figurepath '\RheometricData'],'-depsc2');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[figurepath '\RheometricData'],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
set(gcf, 'color', 'none');
set(findall(gcf,'type','axes'), 'color', 'none');
export_fig RheometricData.png;
movefile('RheometricData.png', figurepath);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');






























