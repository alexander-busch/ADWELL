%% Clean up
clear;
clc;
close all;
global Exp Ref PipeVis CFD Rheo

%% Set path
OriginalDir = pwd;


%% Definitions
DefinePathFiles;
DefineGeometry;


%% Create tables
CreateExpTable;
fluidlist = categories(categorical(Exp.Fluid));
CreateRefTable;
CreatePipeVisTable;
% CreateCFDTable;
% CreateRheoTable;


%% Pre-Plot
fig = figure;
Format_AdjustSize;
Format_Colors;



% Sublot definiton
subplotlist = {'p_dpdl_mdot'; ... & Variables for H2O case
    'p_f_Usl'; ...
    'p_f_ReG'; ...   
    'PipeVisco'; ...
    'a_dpdl_mdot'; ...
    'a_f_Usl'; ...
    'a_f_ReG'; ...
    'Rheometry'};


%% Plot

for i = 1:length(subplotlist)
    subplotname = subplotlist{i}; % Get current subplotname
    subplot(2,4,i); % Activate current subplot
    hold on;
    
    for j = 1:length(fluidlist)
        fluidname = fluidlist(j); % Get current fluid name
        color = colorlist{j}; % Get corresponding plot color
        skip = 0;
        % Set filter in order to...
        if strcmp(fluidname,'PAC1.2') || strcmp(fluidname,'PAC1.5') % skip fluids
            skip = 1;
            
        elseif strcmp(fluidname,'PAC2')
            filter_Exp = strcmp(Exp.Fluid,fluidname) & Exp.RPM==0 & datetime(PipeVis.Date)=='29-Jul-2016'; %strcmp(Exp.Date,'29-Jul-16'); % apply specific filter settings
            filter_Ref = strcmp(Ref.Fluid,fluidname) & Ref.RPM==0; % apply general filter settings for computed data
        
        else
            filter_Exp = strcmp(Exp.Fluid,fluidname) & Exp.RPM==0; % apply general filter settings for experimental data
            filter_Ref = strcmp(Ref.Fluid,fluidname) & Ref.RPM==0; % apply general filter settings for computed data
            
        end
        
        if skip == 1 
            clear skip;
        else
            PlotExpData( subplotname, filter_Exp, filter_Ref, color, fluidname );
        end
    end
    
    FormatSubplots ( subplotname );
    
end


% ReG = ;
% figure
% loglog ( logspace(0,5), FanningFrictionFactor_Ref(logspace(0,5), 'H2O', 3),'--','Color','black');
% 
% round(min(Exp.ReG_p),-2)
% round(max(Exp.ReG_p),-4)
% 












% subplot(2,4,3)
% plot(Exp.ReG_p(filter),Ref.f_p(filter),'color','g')


% figure;
% 
% filter = strcmp(Exp.Fluid,'H2O') & Exp.RPM==0 & ;
% scatter(Exp.MassFlowRate(t), Exp.dpdl_p(t));
% plot(Ref.MassFlowRate(t), Ref.dpdl_p(t));
% 
% filter = strcmp(Exp.Fluid,'PAC1') & Exp.RPM==0;
% scatter(Exp.MassFlowRate(t), Exp.dpdl_p(t));
% plot(Ref.MassFlowRate(t), Ref.dpdl_p(t));
% 
% filter = strcmp(Exp.Fluid,'PAC2') & Exp.RPM==0 & strcmp(Exp.Date,'29-Jul-16');
% scatter(Exp.MassFlowRate(t), Exp.dpdl_p(t));
% plot(Ref.MassFlowRate(t), Ref.dpdl_p(t));
% 
% filter = strcmp(Exp.Fluid,'PAC4') & Exp.RPM==0;
% scatter(Exp.MassFlowRate(t), Exp.dpdl_p(t));
% plot(MeanOfElements(Ref.MassFlowRate(t), 2), MeanOfElements(Ref.dpdl_p(t), 2)); % Take means of every second element
