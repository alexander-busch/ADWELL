close all;
clear all;
clc;


% Read data
load('turbulencedata.mat');

%% Plot format cell arrays

% Turbulence models: ID, LineStyle, LineWidth, Marker, MarkerSize
turbmodels = [unique(cases.Models), {'none'; '-'; '--'; ':'; '-.'}, {2; 1; 1; 1; 1}, {'s'; 'o'; 'o'; 'o'; 'o'}, {10; 4; 4; 4; 4}];

% Fluids: ID, Color
fluids = [unique(cases.Fluid), {'green'; 'cyan'; 'blue'; 'red'}];

% Geometry
geo = unique(cases.Geometry);

% Initialize figure handle array
fig_handle = zeros(2*length(geo),1);


%% Plot data

% Loop all geometries
for aa=1:length(geo)

    % Create dimensional plots dpdx = f(mdot) for pipe and annulus
    fig_handle(2*aa-1) = CreateFigure(['Pressure loss vs. mass flow rate (' geo{aa} ')'],'Mass flow rate dm/dt [kg/s]', 'Pressure gradient dp/dx [Pa/m]', {'lin' 'lin'}, 'DINA5');
    
    % Create non-dimensional plots f = f(Re_G) for pipe and annulus
    fig_handle(2*aa) = CreateFigure(['friction factor vs. generalized Reynolds number (' geo{aa} ')'],'Generalized Reynolds number Re_G [-]', 'Fanning friction factor [-]', {'log' 'log'}, 'DINA5');
    
    % Plot friction factor correlations
    Correlations;
    
    % Create table subset for pipe or annulus    
    cases2print = cases(strcmp(cases.Geometry,geo(aa)) & strcmp(cases.u_x,'x'),:);
    cases2print = sortrows(cases2print,'Description');
    
    % Expand legendentries
    length = length(legendentries);
    legendentries = [legendentries; cell(height(cases2print),1)];
    
    % Loop all table rows
    for bb=1:height(cases2print)
        
        % Get index of current cases for respective plot formatting cell arrays
        idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
        idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

        % dpdx vs. mdot
        figure(fig_handle(2*aa-1));
        plot(cases2print.mdot{bb,1},cases2print.dpdx{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % f vs. Re_G
        figure(fig_handle(2*aa));
        plot(cases2print.Re_G{bb,1},cases2print.f{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % Legendentry
        legendentries{length+bb,1} = [cases2print.Fluid{bb} '; ' cases2print.Models{bb}];
        
    end
    
    % Legend dpdx vs. mdot
    figure(fig_handle(2*aa-1));
    legend(legendentries{:,1}, 'location', 'northwest');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
        
    % Legend f vs. Re_G
    figure(fig_handle(2*aa));
    legend(legendentries{:,1}, 'location', 'northwest');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
end
