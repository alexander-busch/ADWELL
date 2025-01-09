close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
addpath([pwd '\utilities']);
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Turbulence\Figures\';

% Read data
%load('Copy_of_turbulencedata.mat');
 load('turbulencedata.mat');

%% Plot format cell arrays

% Turbulence models: ID, LineStyle, LineWidth, Marker, MarkerSize
turbmodels = [unique(cases.Models), {'none'; '-'; '--'; ':'; '-.'}, {1; 1; 1; 1; 1}, {'s'; 'o'; 'o'; 'o'; 'o'}, {10; 2; 2; 2; 2}];

% Fluids: ID, Color
fluids = sortrows([unique(cases.Fluid), {'green'; 'cyan'; 'blue'; 'red'}, {4; 1; 2; 3}], 3);

% Geometry
geo = unique(cases.Geometry);

% Initialize figure handle array (Rheology, Dim. & non-dim. for pipe and
% annulus)
fig_handle = zeros(1+2*length(geo),1);

% Create sorting index of fluids for plot purposes in case table
looplength = size(fluids); looplength = looplength(1);
for aa=1:height(cases)
    for bb=1:looplength
        if strcmp(cases.Fluid(aa),fluids(bb,1))
            cases.Fluidsortindex{aa} = (fluids{bb,3});
        end
    end
end
% Move variable after fluids, does not work in this release, requires
% presumably R2018
% cases = movevars(cases,'Fluidsortindex','After','Fluid');


%% Plot rheology

% Create rheology table by sorting for fluids and extracting the first
% occurences of each
cases = sortrows(cases,'Fluidsortindex');
[~,out_first,~] = unique(cell2mat(cases.Fluidsortindex), 'first');
rheology = array2table(cases{out_first,{'n_PL_scaled','K_PL_scaled','n_Cr_scaled','K_Cr_scaled','eta_0_scaled','eta_inf_scaled'}});
rheology.Fluid = unique(cases.Fluid,'stable');
rheology.Properties.VariableNames = {'n_PL_scaled','K_PL_scaled','n_Cr_scaled','K_Cr_scaled','eta_0_scaled','eta_inf_scaled','Fluid'};

clear out_first;

% Create rheology plot
fig_handle(1) = CreateFigure('Viscosity vs. shear rate for fluids used','Shear rate d\gamma/dt [1/s]', 'Viscosity \eta [Pa.s]', {'log' 'log'}, 'DINA5');

% Shear rate
gammadot = logspace(-2,4)';

% Loop all table  rows (=fluids)
for aa=1:height(rheology)
    
    % Compute viscosities  
    if strcmp(rheology.Fluid{aa},'Cross')
       rheology.data{aa,1} = [gammadot, rheology.eta_inf_scaled(aa)+(rheology.eta_0_scaled(aa)-rheology.eta_inf_scaled(aa))./(1+(rheology.K_Cr_scaled(aa).*gammadot).^(1-rheology.n_Cr_scaled(aa)))];
    else
       rheology.data{aa,1} = [gammadot, rheology.K_PL_scaled(aa).*gammadot.^(rheology.n_PL_scaled(aa)-1)];
    end

    % Limiting case of Cross model for mu_inf = 0 and large shear rates, i.e. (lambda_Cr gamma_dot)^(1-n_Cr) >> 1
    % K = rheology.eta_0_scaled(aa)*rheology.K_Cr_scaled(aa)^(rheology.n_Cr_scaled(aa)-1);
    % plot(gammadot,K.*gammadot.^(rheology.n_Cr_scaled(aa)-1),'-k');
    
    % Get index of current fluid for from formatting cell arrays
    idx_fluid = find(strcmp(fluids,rheology.Fluid(aa)));
    
    plot(rheology.data{aa,1}(:,1), rheology.data{aa,1}(:,2),...
        'LineStyle','-',...
        'Color',fluids{idx_fluid,2},...
        'LineWidth',2);
    
end

clear gammadot;

% Legend
legend(fluids{:,1}, 'location', 'northeast');
    
% Adjust axis
set(gca, 'ylim', [1e-5 1e-3]);

TightFigure(1);
        
% Print to files
fig_name = 'FlowCurves';
set(fig_handle(1),'PaperPositionMode','auto') %set paper pos for printing
print(fig_handle(1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'



%% Plot data

% Loop all geometries
looplength = size(geo); looplength = looplength(1);
for aa=1:looplength
    % aa=2; % Debugging
    
    % Create dimensional plots dpdx = f(mdot) for pipe and annulus
    fig_handle(1+2*aa-1) = CreateFigure(['Pressure gradient vs. mass flow rate (' geo{aa} ')'],'Mass flow rate dm/dt [kg/s]', 'Pressure gradient dp/dx [Pa/m]', {'lin' 'lin'}, 'DINA5');
    
    % Create non-dimensional plots f = f(Re_G) for pipe and annulus
    fig_handle(1+2*aa) = CreateFigure(['friction factor vs. generalized Reynolds number (' geo{aa} ')'],'Generalized Reynolds number Re_G [-]', 'Fanning friction factor f [-]', {'log' 'log'}, 'DINA5'); % 'PowerPoint'
    
    % Plot friction factor correlations
    Correlations;
    
    % Create table subset for pipe or annulus    
    cases2print = cases(strcmp(cases.Geometry,geo(aa)) & strcmp(cases.u_x_spec,'x'),:);
    
    if strcmp(geo(aa),'Annulus')
        %cases2print = cases2print(cell2mat(cases2print.e)==0 & contains(cases2print.Description,'rpm=0')==1,:);
        cases2print = cases2print(cases2print.e==0 & contains(cases2print.Description,'rpm=0')==1,:);
    end
    
    cases2print = sortrows(cases2print,{'Fluidsortindex','Turbulence','LowRe','nNlowReMod'});
    
    
    
    % Expand legendentries
    legendentrieslength = size(legendentries); legendentrieslength = legendentrieslength(1);
    legendentries = [legendentries; cell(height(cases2print),1)];
    
    % Loop all table rows
    for bb=1:height(cases2print)
        % bb=10;
        % Get index of current cases for respective plot formatting cell arrays
        idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
        idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

        % dpdx vs. mdot
        figure(fig_handle(1+2*aa-1));
        plot(cases2print.mdot{bb,1},cases2print.dpdx{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % f vs. Re_G
        figure(fig_handle(1+2*aa));
        plot(cases2print.Re_G{bb,1},cases2print.f{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % figure
        % plot(cases2print.Re_G{bb,1},cases2print.f{bb,1},':k');
            
        
        
        
        % Legendentry
        legendentries{legendentrieslength+bb,1} = [cases2print.Fluid{bb} '; ' cases2print.Models{bb}];
        
    end
    
    % Legend dpdx vs. mdot
    figure(fig_handle(1+2*aa-1));
    legend(legendentries{5:end,1}, 'location', 'northwest');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);
        
    % Legend f vs. Re_G
    figure(fig_handle(1+2*aa));
    legend(legendentries{:,1}, 'location', 'northeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);
    
    
    % Print to files
    fig_name = [geo{aa} '_dpdx-mdot'];
    set(fig_handle(1+2*aa-1),'PaperPositionMode','auto') %set paper pos for printing
    print(fig_handle(1+2*aa-1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(1+2*aa-1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(1+2*aa-1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
            
    fig_name = [geo{aa} '_f-Re'];
    set(fig_handle(1+2*aa),'PaperPositionMode','auto') %set paper pos for printing
    print(fig_handle(1+2*aa),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(1+2*aa),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(1+2*aa),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'

end


%% Create combined u+y+-plot

% fig_handle(length(fig_handle)) = CreateFigure('u+ vs. y+ for dp/dx = 0.01447 [Pa/m]','y+ [-]', 'u+ [-]', {'log' 'lin'}, 'DINA5');
% 
% % Initialize legendentries
% legendentries = cell(height(cases2print),1);
% 
% for aa = 1:height(cases2print)
%     
%     %Skip case with weird velocity profile
%     if aa~=3
% 
%         bb=6; % dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04];
% 
%         % Get index of current cases for respective plot formatting cell arrays
%         idx_turb = find(strcmp(turbmodels,cases2print.Models(aa)));
%         idx_fluid = find(strcmp(fluids,cases2print.Fluid(aa)));
% 
% 
%         plot(cases2print.y_plus{aa,1}(bb),cases2print.u_plus{aa,1}(bb),...
%             'LineStyle',turbmodels{idx_turb,2},...
%             'Color',fluids{idx_fluid,2},...
%             'LineWidth',turbmodels{idx_turb,3},...
%             'Marker',turbmodels{idx_turb,4},...
%             'MarkerSize',turbmodels{idx_turb,5},...
%             'MarkerEdgeColor',fluids{idx_fluid,2},...
%             'MarkerFaceColor','none');
% 
%         % Legendentry
%         legendentries{aa,1} = [ cases2print.Fluid{aa} '; ' cases2print.Models{aa} ];
%     end
% end 
% 
% % Plot Newtonian theory
% y_plus = linspace(0,12);
% plot(y_plus,y_plus,'-k');
% text(max(y_plus),max(y_plus), ' Viscous sublayer (u+ = y+)','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10,'FontWeight','normal','Color','k','Rotation',65); %'interpreter','latex'
% 
% y_plus = linspace(9,max(xlim)); u_plus = log(y_plus)./0.41+5;
% plot(y_plus,u_plus,'-k');
% text(min(y_plus),min(u_plus), 'Log law (u+ = ln(y+)/0.41+5) ','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10,'FontWeight','normal','Color','k','Rotation',20); %'interpreter','latex'
% 
% % Legend
% legendentries(3) = [];
% legend(legendentries{:,1}, 'location', 'southeast');
% % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
% 
% % Format figure
% set(gca, 'xlim', [1 300],...
%     'XScale', 'log',...  
%     'YScale', 'lin');
% TightFigure(1);




%% Create individual velocity profile and u+ y+ plots

% Get length of figure handle in order to add to it in loop
length_fig_handle = length(fig_handle);

% Loop all fluids
looplength = size(fluids); looplength = looplength(1);
for aa=1:looplength
    % aa=3; % Debugging

    % Create dimensional plots f = f(Re_G)
    fig_handle(length_fig_handle+2*aa-1) = CreateFigure(['Velocity profile for ' fluids{aa} ' and dp/dx = 0.01447 [Pa/m]'],'y/R [-]', 'u(y)/U [-]', {'lin' 'lin'}, 'DINA5');
  
    % Create non-dimensional u+ y+ plot
    fig_handle(length_fig_handle+2*aa) = CreateFigure(['u+ vs. y+ for ' fluids{aa} ' and dp/dx = 0.01447 [Pa/m]'],'y+ [-]', 'u+ [-]', {'log' 'lin'}, 'DINA5');
    
    
    
    % Create table subset for current fluid
    cases2print = cases(strcmp(cases.Geometry,'Pipe') & strcmp(cases.Fluid,fluids(aa)) & strcmp(cases.u_x_spec,'x'),:);
    cases2print = sortrows(cases2print,{'Fluidsortindex','Turbulence','LowRe','nNlowReMod'});

    % Initialize legendentries
    legendentries = cell(height(cases2print),1);

    for bb = 1:height(cases2print)
        % bb=3; % Debugging

        %Skip case with weird velocity profile
        if strcmp(cases2print.ID(bb),'180521_1') %     bb~=3
            % Skip 
        else
            cc=6; % dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04];

            % Get index of current cases for respective plot formatting cell arrays
            idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
            idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));
            
            
            % Create dimensional u_x/U = f(y/R) plots
            figure(fig_handle(length_fig_handle+2*aa-1));
            plot(cases2print.y{bb,1}(cc,:)/(cases2print.d_o(bb)/2),cases2print.u_x{bb,1}(cc,:)/cases2print.U_x{bb}(cc),...
                'LineStyle',turbmodels{idx_turb,2},...
                'Color',fluids{idx_fluid,2},...
                'LineWidth',turbmodels{idx_turb,3},...
                'Marker',turbmodels{idx_turb,4},...
                'MarkerSize',turbmodels{idx_turb,5},...
                'MarkerEdgeColor',fluids{idx_fluid,2},...
                'MarkerFaceColor','none');   
            
            
            % Create non-dimensional u+ y+ plot
            figure(fig_handle(length_fig_handle+2*aa));
            plot(cases2print.y_plus{bb,1}(cc,:),cases2print.u_plus{bb,1}(cc,:),...
                'LineStyle',turbmodels{idx_turb,2},...
                'Color',fluids{idx_fluid,2},...
                'LineWidth',turbmodels{idx_turb,3},...
                'Marker',turbmodels{idx_turb,4},...
                'MarkerSize',turbmodels{idx_turb,5},...
                'MarkerEdgeColor',fluids{idx_fluid,2},...
                'MarkerFaceColor','none');

            % Legendentry
            legendentries{bb,1} = [ cases2print.Fluid{bb} '; ' cases2print.Models{bb} ];
        end
    end 

    % Plot Newtonian theory
    y_plus = linspace(0,12);
    plot(y_plus,y_plus,'-k');
    text(max(y_plus),max(y_plus), ' Viscous sublayer (u+ = y+)','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10,'FontWeight','normal','Color','k','Rotation',65); %'interpreter','latex'

    y_plus = linspace(9,max(xlim)); u_plus = log(y_plus)./0.41+5;
    plot(y_plus,u_plus,'-k');
    text(min(y_plus),min(u_plus), 'Log law (u+ = ln(y+)/0.41+5) ','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10,'FontWeight','normal','Color','k','Rotation',20); %'interpreter','latex'

    % Legend
    if aa==1
        legendentries(3) = [];
    end
    legend(legendentries{:,1}, 'location', 'southeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])

    % Format figure
    set(gca, 'xlim', [1 300],...
        'XScale', 'log',...  
        'YScale', 'lin');
    TightFigure(1);
    
    figure(fig_handle(length_fig_handle+2*aa-1));
    TightFigure(1);
      
    % figure(fig_handle(length_fig_handle+2*aa-1));
    axishandle = findobj(fig_handle(length_fig_handle+2*aa),'type','axes');
    copiedaxishandle = copyobj(axishandle,findobj(fig_handle(length_fig_handle+2*aa-1),'type','figure'));
    set(copiedaxishandle, 'Position',[0.15 0.2 0.75 0.55]);
    legendhandle = axishandle.Legend;
    % legend(legendhandle.String);
    legend(copiedaxishandle,legendhandle.String,'Location','Southeast'); % 'Position',[0.62 0.22 legendhandle.Position(3) legendhandle.Position(4)])
    clear axishandle legendhandle;

        
	% Print to files
    geo = 'Pipe'; % Loop geo if annulus data is available
    fig_name = [geo '_VelocityProfiles_' fluids{aa}];
    set(fig_handle(length_fig_handle+2*aa-1),'PaperPositionMode','auto') %set paper pos for printing
    print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'

end


%% Create dimensional plots U = f(dpdx) for pipe and annulus
fig_handle(length(fig_handle)+1) = CreateFigure(['Pressure gradient sweep for e = 0.5 and 0 rpm'],'Pressure gradient [Pa/m]', 'Mass flow rate dm/dt [kg/s]', {'lin' 'lin'}, 'DINA5');
    
% Create table subset for PL cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0.5 ...
    & contains(cases.Description,'rpm=0')==1 ...
    & contains(cases.Fluid,'PL')==1 ...
        ,:);
        
% Loop all table rows
for bb=1:height(cases2print)

    % Get index of current cases for respective plot formatting cell arrays
    idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
    idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

    % dpdx vs. mdot
    figure(fig_handle(end));
    plot(cases2print.mdot{bb,1},cases2print.dpdx{bb,1},...
        'LineStyle',turbmodels{idx_turb,2},...
        'Color',fluids{idx_fluid,2},...
        'LineWidth',turbmodels{idx_turb,3},...
        'Marker',turbmodels{idx_turb,4},...
        'MarkerSize',turbmodels{idx_turb,5},...
        'MarkerEdgeColor',fluids{idx_fluid,2},...
        'MarkerFaceColor','none');
end

% Create table subset for Newtonian cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0.5 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);
    
    cases2print.rpm{2}(1)
    
    
%% Create dimensional plots U = f(dpdx) for annulus

eccentricity = [0 0.25 0.5 0.75];

for aa = 1:4
    % aa = 2;
    
    
    
% Create table subset for PL cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==eccentricity(aa) ...
    & contains(cases.Description,'rpm=sweep')==1 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);

fig_handle(length(fig_handle)+1) = CreateFigure(['Pressure gradient sweep for e = 0.5 and 0 rpm'],'Pressure gradient [Pa/m]', 'Mass flow rate dm/dt [kg/s]', {'lin' 'lin'}, 'DINA5');

    
% Loop all table rows
for bb=1:height(cases2print)

    % Get index of current cases for respective plot formatting cell arrays
    idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
    idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

    % dpdx vs. mdot
    figure(fig_handle(end));
    plot(cases2print.rpm{bb,1},cases2print.U_x{bb,1},...
        'LineStyle',turbmodels{idx_turb,2},...
        'Color',fluids{idx_fluid,2},...
        'LineWidth',turbmodels{idx_turb,3},...
        'Marker',turbmodels{idx_turb,4},...
        'MarkerSize',turbmodels{idx_turb,5},...
        'MarkerEdgeColor',fluids{idx_fluid,2},...
        'MarkerFaceColor','none');
end

% Create table subset for Newtonian cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0.5 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);
    
    cases2print.rpm{2}(1)
    
    
    