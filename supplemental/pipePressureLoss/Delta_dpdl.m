figure;
hold on;

for j = 1:length(fluidlist)
        fluidname = fluidlist(j); % Get current fluid name
        color = colorlist{j}; % Get corresponding plot color
        skip = 0;
        % Set filter in order to...
        if strcmp(fluidname,'PAC1.2') || strcmp(fluidname,'PAC1.5') % skip fluids
            skip = 1;
            
        elseif strcmp(fluidname,'PAC2')
            filter = strcmp(Exp.Fluid,fluidname) & Exp.RPM==0 & strcmp(Exp.Date,'29-Jul-16'); % apply specific filter settings
        
        else
            filter = strcmp(Exp.Fluid,fluidname) & Exp.RPM==0; % apply general filter settings
            
        end
        
        if skip == 1 
            clear skip;
        else
            % scatter(Exp.MassFlowRate(filter),Exp.delta_dpdl_p(filter),'MarkerEdgeColor',color);
            % scatter(Exp.MassFlowRate(filter),Exp.delta_dpdl_p_normalized(filter),'MarkerEdgeColor',color);
            scatter(Exp.MassFlowRate(filter),Exp.dpdl_p_ratio(filter),'MarkerEdgeColor',color);
        end
end


%% Format

set(gcf,...
    'color','w');

grid('on');

set(gca,...
    'box','on',...
    'FontSize',18,...
    'xlim', [0 1.2],...
    'ylim', [1 2]);

xlabel('Mass flow rate [kg/s]');
ylabel('Delta dpdl_p = dpdl_p_-_E_x_p - dpdl_p_-_C_o_m_p  [Pa/m]');
ylabel('Delta dpdl_p normalized with dpdl_p_-_E_x_p [-]');
ylabel('Delta dpdl_p normalized with dpdl_p_-_C_o_m_p [-]');
ylabel('Delta dpdl_p / dpdl_p_-_C_o_m_p [-]');
  legend('H2O', 'PAC1', 'PAC2', 'PAC4');
    
%  
%     'XScale','log',...
%     'YScale','log',...
%     'xlim', [10 100000],...
%     'ylim', [0.001 1],...
%     'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
%     'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...

