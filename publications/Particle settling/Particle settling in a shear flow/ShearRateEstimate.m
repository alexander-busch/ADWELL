% Comparison of shear rates

% Create shear rate matrix
SR = zeros(length(VolFlowRate), length(d_h));

% Loop all hydraulic diameters
for i = 1:length(d_h)
    
    % Create legend cell array
    Legend = cell(2*length(fluidlist),1);
    
    % Shear rate at the wall estimate
    U_f = VolFlowRate./A(i);
    SR_mean = 12.*U_f./d_h(i);
    SR(i,:) = SR_mean';
    
    CreateFigure(['Shear rate comparison for d_h = ' num2str(d_h(i),2) ' m'],...
        'Reynolds number Re_d_h [-]','Shear rate [1/s]');
    
    % Loop all fluids
    for j = 1:length(fluidlist)
        
        % Fluid apparent viscosity based on wall shear rate estimate
        eta = CrossCoefficients.mu_inf{j}+(CrossCoefficients.mu_0{j}...
            -CrossCoefficients.mu_inf{j})./...
            (1+(CrossCoefficients.lambda{j}.*SR_mean).^CrossCoefficients.n{j});

        % Reynolds number
        Re = rho_f.*U_f.*d_h(i)./eta;
        
        % Fluctuating contribution to total shear rate
        SR_fluc = 0.01107539166.*SR_mean.*Re.^(5/16); % Acc. to Fluent R17.2 User Guide, page 287ff
        % SR_fluc = 0.08333333333.*SR_f.*Re.^(1/2); % Acc. to CFD online, https://www.cfd-online.com/Wiki/Introduction_to_turbulence/Turbulence_kinetic_energy
                
        % Plot
        %scatter(SR_f(i),SR_tur,'d','MarkerEdgeColor',colorlist{j});
        plot(Re,SR_mean,'-d',...
            'Color',colorlist{j},...
            'MarkerEdgeColor',colorlist{j});
        Legend{2*j-1} = ['SR_m_e_a_n (' fluidlist{j} ')'];
        
        plot(Re,SR_fluc,'--^',...
            'Color',colorlist{j},...
            'MarkerEdgeColor',colorlist{j});
        Legend{2*j} = ['SR_f_l_u_c (' fluidlist{j} ')'];
    end

    % Display legend
    Legend = legend(Legend);
    set(Legend,'Location','northwest');
    
end
    