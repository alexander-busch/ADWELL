close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
addpath([pwd '\utilities']);
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Turbulence\Figures\';

% Read data
 load('turbulencedata.mat');
 

%% Plot format cell arrays


% Turbulence models: ID, LineStyle, LineWidth, Marker, MarkerSize
turbmodels = [unique(cases.Models), {'none'; 'none'; '-'; '--'; ':'; '-.'}, {1; 1; 1; 1; 1; 1}, {'s'; 's'; 'o'; 'o'; 'o'; 'o'}, {10; 10; 2; 2; 2; 2}];

% Fluids: ID, Color
fluids = sortrows([unique(cases.Fluid), {'green'; 'cyan'; 'blue'; 'red'}, {4; 1; 2; 3}], 3);

% Geometry
geo = unique(cases.Geometry); geo = flip(geo);

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


%% Exclude data from all plots

% Remove k-omega SST low-Re
toDelete = (strcmp(cases.Turbulence,'k-omega SST') & strcmp(cases.LowRe,'Low Re'));
cases(toDelete,:) = [];
turbmodels = turbmodels(1:5,:);

% Remove n=0.75 pipe case 
toDelete = (strcmp(cases.Geometry,'Pipe') & strcmp(cases.Fluid,'PL') & cases.n_PL_scaled==0.75);
cases(toDelete,:) = [];


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
legend({'H_2O', 'Newtonian #2', 'PL #1', 'Cross #1'}, 'location', 'northeast');
    
% Adjust axis
set(gca, 'ylim', [1e-5 1e-3]);

TightFigure(1);
        
% Print to files
fig_name = 'FlowCurves';
print(fig_handle(1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig_handle(1),[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Plot data

% Loop all geometries
looplength = size(geo); looplength = looplength(1);
for aa=1:looplength
    % aa=2; % Debugging
    
    % Create dimensional plots dpdx = f(mdot) for pipe and annulus
    fig_handle(1+2*aa-1) = CreateFigure(['Pressure gradient vs. mass flow rate (' geo{aa} ')'],'Mass flow rate dm/dt [kg/s]', 'Pressure gradient dp/dx [Pa/m]', {'lin' 'lin'}, 'DINA5');
    
    % Create non-dimensional plots f = f(Re_MR) for pipe and annulus
    fig_handle(1+2*aa) = CreateFigure(['friction factor vs. Metzner-Reed Reynolds number (' geo{aa} ')'],'Metzner-Reed Reynolds number Re_M_R [-]', 'Fanning friction factor f [-]', {'log' 'log'}, 'DINA5'); % 'PowerPoint'
    
    % Plot friction factor correlations on f-Re_MR plot
    Correlations;
    
    % Create table subset for pipe or annulus    
    cases2print = cases(strcmp(cases.Geometry,geo(aa)) & strcmp(cases.u_x_spec,'x'),:);
    
    % Further filtering of cases2print
    if strcmp(geo(aa),'Annulus')
        cases2print = cases2print(cases2print.e==0 & contains(cases2print.Description,'rpm=0')==1,:);
    else

    end
    
    % Sort rows
    cases2print = sortrows(cases2print,{'Fluidsortindex','Turbulence','LowRe','nNlowReMod'});
    
    % Expand legendentries (created in "Correlations"
    legendentrieslength = size(legendentries); legendentrieslength = legendentrieslength(1);
    legendentries = [legendentries; cell(height(cases2print),1)];
      
    % Loop all table rows
    for bb=1:height(cases2print)
        % bb=4;
        % Get index of current cases for respective plot formatting cell arrays
        idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
        idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

        if strcmp(cases2print.Models{bb,1},'DNS2')
           
                      
        else
        
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
        
        % f vs. Re_MR
        figure(fig_handle(1+2*aa));
        plot(cases2print.Re_MR{bb,1},cases2print.f{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % figure
        % plot(cases2print.Re_MR{bb,1},cases2print.f{bb,1},':k');
        
        % Legendentry
        legendentries{legendentrieslength+bb,1} = [cases2print.Fluid{bb} '; ' cases2print.Models{bb}];
        end
    end
    
        legendentries = legendentries(~cellfun('isempty',legendentries));
    
    % Legend dpdx vs. mdot
    figure(fig_handle(1+2*aa-1));
    legend(legendentries{5:end,1}, 'location', 'northwest');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);
        
    % Legend f vs. Re_MR
    figure(fig_handle(1+2*aa));
    legend(legendentries{:,1}, 'location', 'northeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);
    
    
    fig_name = [geo{aa} '_dpdx-mdot'];
    print(fig_handle(1+2*aa-1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(1+2*aa-1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(1+2*aa-1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
    print(fig_handle(1+2*aa-1),[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
    
    fig_name = [geo{aa} '_f-Re_MR'];
    print(fig_handle(1+2*aa),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(1+2*aa),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(1+2*aa),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
    print(fig_handle(1+2*aa),[fig_path fig_name],'-dsvg');   
    
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
    % aa=2; % Debugging
    % aa=3; % Debugging

    % Create dimensional plots f = f(Re_MR)
    fig_handle(length_fig_handle+2*aa-1) = CreateFigure(['Velocity profile for ' fluids{aa} ' and dp/dx = 0.01447 [Pa/m]'],'y/R [-]', 'u(y)/U [-]', {'lin' 'lin'}, 'DINA5');
  
    % Create non-dimensional u+ y+ plot
    fig_handle(length_fig_handle+2*aa) = CreateFigure(['u+ vs. y+ for ' fluids{aa} ' and dp/dx = 0.01447 [Pa/m]'],'y+ [-]', 'u+ [-]', {'log' 'lin'}, 'DINA5');
    
    
    
    % Create table subset for current fluid
    cases2print = cases(strcmp(cases.Geometry,'Pipe') & strcmp(cases.Fluid,fluids(aa)) & strcmp(cases.u_x_spec,'x'),:);
    
    % Sort
    cases2print = sortrows(cases2print,{'Fluidsortindex','Turbulence','LowRe','nNlowReMod'});

    % Initialize legendentries and max yplus value
    legendentries = cell(height(cases2print),1);
    yplus_max = 0;
    
    for bb = 1:height(cases2print)
        % bb=3; % Debugging

        %Skip case with weird velocity profile
        if (strcmp(cases2print.ID(bb),'180521_1') || contains(cases2print.Models{bb,1},'DNS2')==1)%     bb~=3
            % Skip 
        else
            
            % Design point dpdx = 0.01447 Pa/m
            cc=find(cases2print.dpdx{bb}==0.01447); % dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04];
            
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
            
            yplus_max = max(yplus_max, max((cases2print.y_plus{bb,1}(cc,:))));
            
            yplus = cases2print.y_plus{bb,1}(cc,:);
            yplus = (yplus(~isnan(yplus)));
            uplus = cases2print.u_plus{bb,1}(cc,:);
            uplus = (uplus(~isnan(uplus)));
            uplus_int(bb) = 2*trapz(yplus,uplus.*yplus) / max(yplus).^2 *cases2print.ufricf{bb,1}(cc);
           % text(max(yplus),max(uplus), num2str(uplus_int(bb),3),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontWeight','normal','Color','k','Rotation',0); %'interpreter','latex'

            % Legendentry
            legendentries{bb,1} = [ cases2print.Fluid{bb} '; ' cases2print.Models{bb}  '; Integral = ' num2str(uplus_int(bb),3) ];
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
    if (aa==1 && length(legendentries)>2)
        legendentries(3) = [];
    end
    
    legendentries = legendentries(~cellfun('isempty',legendentries));
    
    legend(legendentries{:,1}, 'location', 'southeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])

    % Format figure
    set(gca, 'xlim', [1 yplus_max],...
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
    fig_name = [geo{1} '_VelocityProfiles_' fluids{aa}];
    print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
    print(fig_handle(length_fig_handle+2*aa-1),[fig_path fig_name],'-dsvg');
end


%% Create dimensional plot dpdx = f(mdot) for Newtonian flow in concentric annulus

cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0 ...
    & contains(cases.Description,'rpm=0')==1 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);

fig_handle(length(fig_handle)+1) = CreateFigure(['dpdx = f(mdot) for Newtonian flow in concentric annulus'], 'Mass flow rate dm/dt [kg/s]', 'Pressure gradient [Pa/m]', {'lin' 'lin'}, 'DINA5');

% Change definition of plot formating,
% such that line is displayed for DNS % data
% turbmodels{1,2} = '-';

% Loop all table rows
for bb=1:height(cases2print)
    % bb=3

    % Get index of current cases for respective plot formatting cell arrays
    idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
    idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

    % mdot vs. rpm
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

% Reset definition of plot formating
%turbmodels{1,2} = 'none';




% Newtonian correlations
bb = 2;

% Blasius
linetype = '--';
mdot = xlim; mdot = (mdot(1):0.5:mdot(2))';
vel = mdot./cases2print.rho_scaled(bb)./cases2print.A{bb};
Re = cases2print.rho_scaled(bb).*cases2print.d_h(bb).*vel./cases2print.K_PL_scaled(bb);
f=0.316.*Re.^-0.25/4;
dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
plot(mdot,dpdx,'LineStyle',linetype,'Color',fluids{idx_fluid,2});
%legendentries{5,1} = 'Blasius w/ Haciislamoglu & Cartalos (1994) correction';

% % Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% f_Morrison = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
% dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
% dpdx = EccentricPressureCorrection( dpdx, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_scaled(aa)  );
% plot(vel,dpdx,'LineStyle',linetype);
% legendentries{7,1} = [ 'Morrison, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm(aa)) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Haaland
linetype = '-.';
f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(bb)).^1.11+(6.9./Re))).^2);
dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
plot(mdot,dpdx,'LineStyle',linetype,'Color',fluids{idx_fluid,2});
%legendentries{6,1} = 'Haaland w/ Haciislamoglu & Cartalos (1994) correction';

TightFigure(1);





%% Create dimensional plot mdot = f(rpm) for fixed dpdx for annulus

cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & contains(cases.Description,'rpm=sweep')==1 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);

fig_handle(length(fig_handle)+1) = CreateFigure(['rpm sweep for fixed dpdx'],'Inner pipe rotation rate [rpm]', 'Mass flow rate dm/dt [kg/s]', {'lin' 'lin'}, 'DINA5');

% Change definition of plot formating,
% such that line is displayed for DNS % data
turbmodels{1,2} = '-';


% Initialize legendentries and max yplus value
legendentries = cell(height(cases2print),1);

% Define colorband

ColorSet = ColorBand(height(cases2print)/2);


% Loop all table rows
for bb=1:height(cases2print)
    
    
    % reset color after fourth plot
    if bb==5
        set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-4);
    end

    % Get index of current cases for respective plot formatting cell arrays
    idx_turb = find(strcmp(turbmodels,cases2print.Models(bb)));
    idx_fluid = find(strcmp(fluids,cases2print.Fluid(bb)));

    % mdot vs. rpm
    figure(fig_handle(end));
    lineH = plot(cases2print.rpm{bb,1},cases2print.mdot{bb,1},...
        'LineStyle',turbmodels{idx_turb,2},...
        'LineWidth',turbmodels{idx_turb,3},...
        'Marker',turbmodels{idx_turb,4},...
        'MarkerSize',turbmodels{idx_turb,5},...
        'MarkerFaceColor','none');
%      'MarkerEdgeColor',fluids{idx_fluid,2},...
%          'Color',fluids{idx_fluid,2},...

    
        % Add e as text to data
    if contains(cases2print.Models{bb,1},'DNS')==1
        if isnan(cases2print.mdot{bb,1})==0
            idx = length(cases2print.mdot{bb,1});
        else           
            idx = find(isnan(cases2print.mdot{bb,1}),1)-1;
        end
        text(...
            cases2print.rpm{bb,1}(idx),...
            cases2print.mdot{bb,1}(idx),...
            ['  ' num2str(round(cases2print.e(bb),2))],...
            'HorizontalAlignment','left',...
            'VerticalAlignment','middle',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',get(lineH, 'Color')); %'interpreter','latex'
        % fluids{idx_fluid,2}
    else
        if isnan(cases2print.mdot{bb,1})==0
            idx = length(cases2print.mdot{bb,1});
        else           
            idx = find(isnan(cases2print.mdot{bb,1}),1)-1;
        end
        
        text(...
            cases2print.rpm{bb,1}(idx),...
            cases2print.mdot{bb,1}(idx),...
            ['  ' num2str(round(cases2print.e(bb),2))],...
            'HorizontalAlignment','left',...
            'VerticalAlignment','middle',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',get(lineH, 'Color')); %'interpreter','latex'
            % fluids{idx_fluid,2}
    end

    
    if (bb==1 || bb==5)
        % Legendentry
        legendentries{bb,1} = [ cases2print.Fluid{bb} '; ' cases2print.Models{bb} ];
    else
        set(get(get( lineH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end    
    
end

% Reset definition of plot formating
turbmodels{1,2} = 'none';

TightFigure(1);


% legendentries = legendentries([1 5],1)
% legend(legendentries{:,1}, 'location', 'southwest');
   
fig_name = [geo{2} '_mdot-rpm'];
print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(length(fig_handle)),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dsvg');

%% Create dimensional plot dpdx = f(mdot) for PL and e = 0.5 for annulus

    
% Create table subset for PL cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0.5 ...
    & contains(cases.Fluid,'PL')==1 ...
        ,:);
    
% Remove k-epsilon Low Re nN mod
toDelete = (strcmp(cases.Turbulence,'k-epsilon Standard') & strcmp(cases.LowRe,'Low Re') & strcmp(cases.LowRe,'nN mod'));
cases(toDelete,:) = [];
turbmodels = turbmodels(1:5,:);
    
fig_handle(length(fig_handle)+1) = CreateFigure(['Pressure gradient sweep for e = 0.5'], 'Mass flow rate dm/dt [kg/s]', 'Pressure gradient [Pa/m]', {'lin' 'lin'}, 'DINA5');
 
% Change definition of plot formating,
% such that line is displayed for DNS % data
% turbmodels{1,2} = '-';
        
% Loop all table rows
for bb=1:8

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
    
    % Add rpm as text to data
    if contains(cases2print.Models{bb,1},'DNS')==1
        text(...
            cases2print.mdot{bb,1}(4),...
            cases2print.dpdx{bb,1}(4),...
            num2str(round(cases2print.rpm{bb,1},1)),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',fluids{idx_fluid,2}); %'interpreter','latex'
    else
        text(...
            cases2print.mdot{bb,1}(6),...
            cases2print.dpdx{bb,1}(6),...
            num2str(round(cases2print.rpm{bb,1},1)),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',fluids{idx_fluid,2}); %'interpreter','latex'
    end
end


% Reset definition of plot formating
% turbmodels{1,2} = 'none';


% Friction factor correlations
bb = 4;

mdot = xlim; mdot = (mdot(1):0.5:mdot(2))';
vel = mdot./cases2print.rho_scaled(bb)./cases2print.A{bb};
Re = cases2print.rho_scaled(bb).*vel.^(2-cases2print.n_PL_scaled(bb)).*cases2print.d_h(bb).^cases2print.n_PL_scaled(bb)./(12.^(cases2print.n_PL_scaled(bb)-1).*cases2print.K_PL_scaled(bb).*((2.*cases2print.n_PL_scaled(bb)+1)./(3.*cases2print.n_PL_scaled(bb))).^cases2print.n_PL_scaled(bb));
Re_cr = 2100*(4*cases2print.n_PL_scaled(bb)+2).*(5*cases2print.n_PL_scaled(bb)+3)./(3*(3*cases2print.n_PL_scaled(bb)+1).^2);
vel_cr = (Re_cr.*(12.^(cases2print.n_PL_scaled(bb)-1).*cases2print.K_PL_scaled(bb).*((2.*cases2print.n_PL_scaled(bb)+1)./(3.*cases2print.n_PL_scaled(bb))).^cases2print.n_PL_scaled(bb))./(cases2print.rho_scaled(bb).*cases2print.d_h(bb).^cases2print.n_PL_scaled(bb))).^(1./(2-cases2print.n_PL_scaled(bb)));

% PL, turbulent regime (Irvine 1988)
linetype = '-';
f=zeros(length(Re),1);
for ii = 1:length(Re)
   if Re(ii)<Re_cr
       f(ii)=24/Re(ii);
   else
       f(ii)=((((2.^(cases2print.n_PL_scaled(bb)+4))./(7.^(7.*cases2print.n_PL_scaled(bb)))).*(4.*cases2print.n_PL_scaled(bb)./(3.*cases2print.n_PL_scaled(bb)+1)).^(3.*cases2print.n_PL_scaled(bb).^2))./Re(ii)).^(1./(3.*cases2print.n_PL_scaled(bb)+1));
   end
end
f_cr = 24/Re_cr;
%f_cr = ((((2.^(cases2print.n_PL_scaled+4))./(7.^(7.*cases2print.n_PL_scaled))).*(4.*cases2print.n_PL_scaled(bb)./(3.*cases2print.n_PL_scaled(bb)+1)).^(3.*cases2print.n_PL_scaled(bb).^2))./Re_cr).^(1./(3.*cases2print.n_PL_scaled(bb)+1));
dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
plot(mdot,dpdx,'LineStyle',linetype,'Color','k');
%legendentries{5,1} = ['Irvine (1988), e = ' num2str(cases2print.e) ', r = 0 m' ];

% PL, turbulent regime (Dodge & Metzner 1959), implicit in f hence specify
% existing f >> Re >> vel
linetype = '--'; 
for ii = 1:length(Re)
   if f(ii)>f_cr
       Re(ii)=24/f(ii);
   else
       Re(ii) = (10.^((1./sqrt(f(ii))+0.4./(cases2print.n_PL_scaled(bb).^1.2)).*(cases2print.n_PL_scaled(bb).^0.75)./4))./(f(ii).^((2-cases2print.n_PL_scaled(bb))./2));
   end
end
vel = (Re.*(12.^(cases2print.n_PL_scaled(bb)-1).*cases2print.K_PL_scaled(bb).*((2.*cases2print.n_PL_scaled(bb)+1)./(3.*cases2print.n_PL_scaled(bb))).^cases2print.n_PL_scaled(bb))./(cases2print.rho_scaled(bb).*cases2print.d_h(bb).^cases2print.n_PL_scaled(bb))).^(1./(2-cases2print.n_PL_scaled(bb)));
dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
plot(mdot,dpdx,'LineStyle',linetype,'Color','k');
%legendentries{6,1} = ['Dodge & Metzner (1959), e = ' num2str(cases2print.e) ', r = 0 m' ]


% Create table subset for Newtonian cases to print
cases2print = cases(strcmp(cases.Geometry,'Annulus') ...
    & cases.d_o==4 ...
    & cases.e==0.5 ...
    & contains(cases.Fluid,'Newtonian')==1 ...
        ,:);

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
    
    % Add rpm as text to data
    if contains(cases2print.Models{bb,1},'DNS')==1
        text(...
            cases2print.mdot{bb,1},...
            (repelem(cases2print.dpdx{bb,1},length(cases2print.mdot{bb,1})))',...
            num2str(round(cases2print.rpm{bb,1},1)),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',fluids{idx_fluid,2}); %'interpreter','latex'
    else
        text(...
            cases2print.mdot{bb,1},...
            (repelem(cases2print.dpdx{bb,1},length(cases2print.mdot{bb,1})))',...
            num2str(round(cases2print.rpm{bb,1},1)),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'Rotation',0,...
            'FontSize',10,...
            'FontWeight','normal',...
            'Color',fluids{idx_fluid,2}); %'interpreter','latex'
    end
 end  

 
% Format figure
set(gca, ...
    'xlim', [2.5 15],...
    'ylim', [0 0.02],...
    'XScale', 'lin',...  
    'YScale', 'lin');

TightFigure(1);



% Newtonian correlations
bb = 2;

% Blasius
linetype = '-';
mdot = xlim; mdot = (mdot(1):0.5:mdot(2))';
vel = mdot./cases2print.rho_scaled(bb)./cases2print.A{bb};
Re = cases2print.rho_scaled(bb).*cases2print.d_h(bb).*vel./cases2print.K_PL_scaled(bb);
% f=0.316.*Re.^-0.25/4;
% dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
% dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
% plot(mdot,dpdx,'LineStyle',linetype,'Color',fluids{idx_fluid,2});
%legendentries{5,1} = 'Blasius w/ Haciislamoglu & Cartalos (1994) correction';

% % Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
plot(mdot,dpdx,'LineStyle',linetype,'Color',[0.6 0.6 0.6]);
% legendentries{7,1} = [ 'Morrison, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm(aa)) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Haaland
% linetype = '-.';
% f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(bb)).^1.11+(6.9./Re))).^2);
% dpdx = f.*cases2print.rho_scaled(bb).*vel.^2./(cases2print.d_h(bb)./2);
% dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(bb), cases2print.d_o(bb), cases2print.e(bb), cases2print.n_PL_scaled(bb)  );
% plot(mdot,dpdx,'LineStyle',linetype,'Color',fluids{idx_fluid,2});
%legendentries{6,1} = 'Haaland w/ Haciislamoglu & Cartalos (1994) correction';

fig_name = [geo{2} '_dpdx-mdot_e'];
print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(length(fig_handle)),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dsvg');

print(fig_handle(16),[fig_path fig_name],'-dsvg');
figure(fig_handle(16))


%% Alternative scaling where apparent viscosity is based on tau_wall and dpdx

% Pipe
aa=1;

    % Create non-dimensional plots f = f(Re_G) for annulus
    fig_handle(length(fig_handle)+1) = CreateFigure(['friction factor vs. Generalized Reynolds number (' geo{aa} ')'],'Reynolds number Re_G [-]', 'Fanning friction factor f [-]', {'log' 'log'}, 'DINA5'); % 'PowerPoint'
    % Plot friction factor correlations
    Correlations;
    
    % Create table subset for pipe or annulus    
    cases2print = cases(strcmp(cases.Geometry,geo(aa)) & strcmp(cases.u_x_spec,'x'),:);
    
        % Further filtering of cases2print
    if strcmp(geo(aa),'Annulus')
        cases2print = cases2print(cases2print.e==0 & contains(cases2print.Description,'rpm=0')==1,:);
    else

    end
    
    % Reorder rows for plotting
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
        
        % f vs. Re_MR
        figure(length(fig_handle));
        plot(cases2print.Re_G{bb,1},cases2print.f{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % Legendentry
        legendentries{legendentrieslength+bb,1} = [cases2print.Fluid{bb} '; ' cases2print.Models{bb}];
        
    end

    % Legend f vs. Re_MR
    figure(fig_handle(length(fig_handle)));
    legend(legendentries{:,1}, 'location', 'northeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);

    fig_name = [geo{aa} '_f-Re_G'];
    print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length(fig_handle)),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
    print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dsvg');
    
    
% Pipe
aa=2;

    % Create non-dimensional plots f = f(Re_G) for annulus
    fig_handle(length(fig_handle)+1) = CreateFigure(['friction factor vs. Generalized Reynolds number (' geo{aa} ')'],'Reynolds number Re_G [-]', 'Fanning friction factor f [-]', {'log' 'log'}, 'DINA5'); % 'PowerPoint'
    % Plot friction factor correlations
    Correlations;
    
    % Create table subset for pipe or annulus    
    cases2print = cases(strcmp(cases.Geometry,geo(aa)) & strcmp(cases.u_x_spec,'x'),:);
    
        % Further filtering of cases2print
    if strcmp(geo(aa),'Annulus')
        cases2print = cases2print(cases2print.e==0 & contains(cases2print.Description,'rpm=0')==1,:);
    else

    end
    
    % Reorder rows for plotting
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
        
        % f vs. Re_MR
        figure(length(fig_handle));
        plot(cases2print.Re_G{bb,1},cases2print.f{bb,1},...
            'LineStyle',turbmodels{idx_turb,2},...
            'Color',fluids{idx_fluid,2},...
            'LineWidth',turbmodels{idx_turb,3},...
            'Marker',turbmodels{idx_turb,4},...
            'MarkerSize',turbmodels{idx_turb,5},...
            'MarkerEdgeColor',fluids{idx_fluid,2},...
            'MarkerFaceColor','none');
        
        % Legendentry
        legendentries{legendentrieslength+bb,1} = [cases2print.Fluid{bb} '; ' cases2print.Models{bb}];
        
    end

    % Legend f vs. Re_MR
    figure(fig_handle(length(fig_handle)));
    legend(legendentries{:,1}, 'location', 'northeast');
    % reorderLegend([13,14,15,1,3,8,11,2,7,5,4,12,9,10,6])
    TightFigure(1);

    fig_name = [geo{aa} '_f-Re_G'];
    print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dpdf','-bestcreate');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length(fig_handle)),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
    % print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
    print(fig_handle(length(fig_handle)),[fig_path fig_name],'-dsvg');

