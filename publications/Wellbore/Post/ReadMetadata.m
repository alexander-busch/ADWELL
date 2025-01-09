%% Read metadata

clear all;
close all;
clc;

filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive';
filepath_metadata = [filepath_simdata  '\' 'metadata2.xlsx'];

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata,'Sheet','Cases');
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
opts = setvartype(opts,{'dpdx', 'rpm', 'rpm_whirl', 'mdot_f', 'Q_f', 'U_f_x', 'mdot_s', 'Q_s', 'U_s_x', 'CTR_mdot', 'CTR_Q'},'char');
% disp([opts.VariableNames' opts.VariableTypes']);
cases = readtable(filepath_metadata,opts,'Sheet','Cases');

% Remove first ad second row as these contains the units and description
%cases.Properties.VariableUnits
%cases.Properties.VariableDescriptions
cases([1,2],:) = [];

% Conversion from character strings to numerical arrays
for aa=1:height(cases)
    % aa = 36
    cases.dpdx{aa} = str2num(cases.dpdx{aa});
    cases.dpdx_exp{aa} = str2num(cases.dpdx_exp{aa});
    cases.mdot_f{aa} = str2num(cases.mdot_f{aa});
    cases.Q_f{aa} = str2num(cases.Q_f{aa});
    cases.U_f_x{aa} = str2num(cases.U_f_x{aa});
    cases.mdot_s{aa} = str2num(cases.mdot_s{aa});
    cases.Q_s{aa} = str2num(cases.Q_s{aa});
    cases.U_s_x{aa} = str2num(cases.U_s_x{aa});
	cases.U_f_x_exp{aa} = str2num(cases.U_f_x_exp{aa});
    cases.rpm{aa} = str2num(cases.rpm{aa});
    cases.rpm_whirl{aa} = str2num(cases.rpm_whirl{aa});
    
end

save('wellboredata.mat', 'cases');


%% Read CFD data

% Loop all table rows
for aa=1:height(cases)
    % aa = 7;
    % aa = 24;
    % aa = 36;
    % aa = 47;
    
    if strcmp(cases.Solid{aa},'n/a')
       % Single-phase
       
        if (strcmp(cases.flow_spec{aa},'mdot') || strcmp(cases.flow_spec{aa},'U_f'))
            % Flow spec

            if (length(cases.dpdx{aa}) == 1)
                % Do nothing, solution contained in xlsx file
                
            else
                
                % !!! Not required for the time being
%            
%                 fullpath = [filepath_simdata '\' cases.ID{aa} '\solution\monitors\dpdxper-pr-grad'];
%             
%                 % Check for missing file type and rename file if necessary
%                 if exist(fullpath, 'file') == 2
%                     movefile(fullpath,[fullpath '.txt']);
%                 end
% 
%                 [iterations, dpdx] = importMonitoredQuantities([fullpath '.txt']);
%                 for bb = 1:length(cases.U_f_x{aa})
%                     idx = find(iterations==(bb*iterations(end)/length(cases.U_f_x{aa})));
%                     cases.dpdx{aa,1}(bb) = dpdx(min(idx));
%                 end
%                 % Transpose to column vector
%                 cases.dpdx{aa}=cases.dpdx{aa}';
%                 % Calculate volumetric flow rate and bulk velocity    
%                 cases.Q_f{aa} = cases.U_f_x{aa}.*cases.A(aa);
%                 cases.mdot_f{aa} = cases.Q_f{aa}.*cases.rho_f(aa);
% 
% 
%                 %             % Plot reansient pressure monitor data           
%                 %             figure; hold on;
%                 %             plot(iterations, -dpdx)
%                 %             set(gca,'ylim', [0 3500]); 
                
            end
            
           
       else
           % Pressure spec
           if (length(cases.dpdx{aa}) == 1)
               % Do nothing, solution contained in xlsx file
           else

                % Read mass flow monitor file
                fullpath = [filepath_simdata '\' cases.ID{aa} '\solution\monitors\fluid-x-mass-flow-rate.txt'];
                [iterations, m_dot] = importMonitoredQuantities(fullpath);

                % Get average of last 500 values of respective interval
                for bb = 1:length(cases.dpdx{aa})
                    idx = find(iterations==(bb*iterations(end)/length(cases.dpdx{aa})));
                    cases.mdot_f{aa,1}(bb) = mean(m_dot(min(idx)-500):m_dot(min(idx)));
                end

                % Transpose to column vector
                cases.mdot_f{aa}=cases.mdot_f{aa}';

                % Calculate volumetric flow rate and bulk velocity
                cases.Q_f{aa} = cases.mdot_f{aa}./cases.rho_f(aa);
                cases.U_f_x{aa} = cases.Q_f{aa}./cases.A(aa);

                % Plot transient superficial velocity monitor data           
                % figure; hold on;
                % plot(iterations, m_dot)
                % set(gca,'ylim', [0 3500]);

            end
        end

    else
        % Multi-phase
        
        if (strcmp(cases.flow_spec{aa},'mdot') || strcmp(cases.flow_spec{aa},'U_f'))
            % Flow spec
                
            % !!! Not required for the time being
           
           
        else
           % Pressure spec
           if (length(cases.dpdx{aa}) == 1)
               % Do nothing, solution contained in xlsx file
           else
                if (length(cases.rpm{aa}) == 1)
                    
                    % Single rpm case

                    % Read mass flow monitor file
                    fullpath = [filepath_simdata '\' cases.ID{aa} '\solution\monitors\fluid-x-mass-flow-rate.txt'];
                    [iterations, m_dot] = importMonitoredQuantities(fullpath);

                    % Get average of last 500 values of respective interval
                    for bb = 1:length(cases.dpdx{aa})
                        idx = find(iterations==(bb*iterations(end)/length(cases.dpdx{aa})));
                        cases.mdot_f{aa,1}(bb) = mean(m_dot(min(idx)-500):m_dot(min(idx)));
                    end

                    % Transpose to column vector
                    cases.mdot_f{aa}=cases.mdot_f{aa}';

                    % Calculate volumetric flow rate and bulk velocity
                    cases.Q_f{aa} = cases.mdot_f{aa}./cases.rho_f(aa);
                    cases.U_f_x{aa} = cases.Q_f{aa}./cases.A(aa);

                    % Plot transient superficial velocity monitor data           
                    % figure; hold on;
                    % plot(iterations, m_dot)
                    % set(gca,'ylim', [0 120]);
                    
                else
                    
                    % rpm sweep
                    
                    for bb = 1:length(cases.rpm{aa})
                        % bb = 4
                        
                        % Read mass flow monitor file
                        fullpath = [filepath_simdata '\' cases.ID{aa} '\' num2str(bb) '\solution\monitors\fluid-x-mass-flow-rate.txt'];
                        if exist(fullpath, 'file') == 2
                            [timeseries, m_dot] = importMonitoredQuantities(fullpath);
                        else
                            timeseries = 0;
                            m_dot = 0;
                        end
                        
%                         if (strcmp(cases.ID{aa},'181026_1') || strcmp(cases.ID{aa},'181026_2'))
%                         
%                             % Concatenate data from two data sets 
%                             
%                             % Reduce arrays to 90 s of flowtime
%                             if (strcmp(cases.ID{aa},'181026_1') && (bb > 4))
%                                 idx = find(timeseries==120.0);
%                             else
%                                 idx = find(timeseries==90.0);
%                             end
%                             timeseries = timeseries(1:idx);
%                             m_dot = m_dot(1:idx);
%                             
%                             % Get second set of data
%                             fullpath = [filepath_simdata '\' cases.ID{aa} '\' num2str(bb) '\solution\monitors\fluid-x-mass-flow-rate2.txt'];
%                             [timeseries2, m_dot2] = importMonitoredQuantities(fullpath);
%                             
%                             % concatenate data
%                             timeseries = [timeseries; timeseries(end)+timeseries2];
%                             m_dot = [m_dot; m_dot2];
%                             
%                             clear timeseries2 m_dot2;
%                         
%                         end
                        
                        % Write timeseries to table 
                        if bb==1
                            cases.timeseries{aa,1} = [];
                            cases.mdot_f_series{aa,bb} = [];
                            cases.Q_f_series{aa,bb} = [];
                            cases.U_f_x_series{aa,bb} = [];
                            cases.mdot_s_series{aa,bb} = [];
                            cases.Q_s_series{aa,bb}= [];
                            cases.U_s_x_series{aa,bb} = [];
                            cases.timeseries{aa}(:,1) = timeseries;
                        end
                        
                        % Write mass flow series to table
                        cases.mdot_f_series{aa,bb}(:,1) = m_dot;
                        
                        % Calculate volumetric flow rate and bulk velocity
                        cases.Q_f_series{aa,bb} = cases.mdot_f_series{aa,bb}./cases.rho_f(aa);
                        cases.U_f_x_series{aa,bb} = cases.Q_f_series{aa,bb}./cases.A(aa);
                        
                        if m_dot~=0 
                            % Get average of last 5 s values of respective interval
                            for cc = 1:length(cases.dpdx{aa})
                                %idx = find(timeseries==(cc*timeseries(end)/length(cases.dpdx{aa})));
                                % Previous formulation requires exact match
                                % while  following formulation finds nearest
                                % value
                                [ ~, idx ] = min( abs( timeseries-(cc*timeseries(end)/length(cases.dpdx{aa})) ) );
                                if cc==1
                                    cases.mdot_f{aa,bb}(cc) = mean(m_dot(min(idx)-1000:min(idx)));
                                else
                                    cases.mdot_f{aa,bb}(cc) = mean(m_dot(min(idx)-5000:min(idx)));
                                end
                                cases.Q_f{aa,bb}(cc) = cases.mdot_f{aa,bb}(cc)./cases.rho_f(aa);
                                cases.U_f_x{aa,bb}(cc) = cases.Q_f{aa,bb}(cc)./cases.A(aa);
                            end
                         else
%                                 cases.mdot_f{aa,bb}(cc) = 0;
%                                 cases.Q_f{aa,bb}(cc) = 0;
%                                 cases.U_f_x{aa,bb}(cc) = 0;
                                cases.mdot_f{aa,bb} = 0;
                                cases.Q_f{aa,bb} = 0;
                                cases.U_f_x{aa,bb} = 0;
                        end                        
                        
                        % Read solid mass flow monitor file
                        %aa= 47
                        %bb= 1
                        
                        fullpath = [filepath_simdata '\' cases.ID{aa} '\' num2str(bb) '\solution\monitors\solid-x-mass-flow-rate.txt'];
                        if exist(fullpath, 'file') == 2
                            [timeseries, m_dot] = importMonitoredQuantities(fullpath);
                        else
                            timeseries = 0;
                            m_dot = 0;
                        end
                        
%                         hold on
%                         plot(timeseries,m_dot)

%                         if (strcmp(cases.ID{aa},'181026_1') || strcmp(cases.ID{aa},'181026_2'))
%                         
%                             % Concatenate data from two data sets 
%                             
%                             % Reduce arrays to 90 s of flowtime
%                             if (strcmp(cases.ID{aa},'181026_1') && (bb > 4))
%                                 idx = find(timeseries==120.0);
%                             else
%                                 idx = find(timeseries==90.0);
%                             end
%                             timeseries = timeseries(1:idx);
%                             m_dot = m_dot(1:idx);
%                             
%                             % Get second set of data
%                             fullpath = [filepath_simdata '\' cases.ID{aa} '\' num2str(bb) '\solution\monitors\solid-x-mass-flow-rate2.txt'];
%                             [timeseries2, m_dot2] = importMonitoredQuantities(fullpath);
%                             
%                             % concatenate data
%                             timeseries = [timeseries; timeseries(end)+timeseries2];
%                             m_dot = [m_dot; m_dot2];
%                             
%                             clear timeseries2 m_dot2;
%                         
%                         end
                        
                        % Write mass flow series to table
                        cases.mdot_s_series{aa,bb}(:,1) = m_dot;
                        
                        % plot(cases.timeseries{aa}(:,1),cases.mdot_s_series{aa,bb}(:,1))
                        
                        % Calculate volumetric flow rate and bulk velocity
                        cases.Q_s_series{aa,bb} = cases.mdot_s_series{aa,bb}./cases.rho_s(aa);
                        cases.U_s_x_series{aa,bb} = cases.Q_s_series{aa,bb}./cases.A(aa);
                        
                        %plot(cases.timeseries{aa}(:,1),cases.U_s_x_series{aa,bb}(:,1))
                        
                        if m_dot~=0 
                            % Get average of last 5 s values of respective interval
                            for cc = 1:length(cases.dpdx{aa})
                                % cc = 2
                                % cc = 3
                                % idx = find(timeseries==(cc*timeseries(end)/length(cases.dpdx{aa})));
                                % Previous formulation requires exact match
                                % while  following formulation finds nearest
                                % value
                                [ ~, idx ] = min( abs( timeseries-(cc*timeseries(end)/length(cases.dpdx{aa})) ) );
                                if cc==1
                                    cases.mdot_s{aa,bb}(cc) = mean(m_dot(min(idx)-1000:min(idx)));
                                else
                                    cases.mdot_s{aa,bb}(cc) = mean(m_dot(min(idx)-5000:min(idx)));
                                end
                                cases.Q_s{aa,bb}(cc) = cases.mdot_s{aa,bb}(cc)./cases.rho_s(aa);
                                cases.U_s_x{aa,bb}(cc) = cases.Q_s{aa,bb}(cc)./cases.A(aa);
                            end
                        else
%                                 cases.mdot_s{aa,bb}(cc) = 0;
%                                 cases.Q_s{aa,bb}(cc) = 0;
%                                 cases.U_s_x{aa,bb}(cc) = 0;

                                cases.mdot_s{aa,bb} = 0;
                                cases.Q_s{aa,bb} = 0;
                                cases.U_s_x{aa,bb} = 0;
                        end
                        
                    end

                    % Compute CTR
                    for bb = 1:length(cases.rpm{aa})
                            cases.CTR_mdot{aa,bb} = cases.mdot_s{aa,bb}./cases.mdot_f{aa,bb};
                            cases.CTR_Q{aa,bb} = cases.U_s_x{aa,bb}./cases.U_f_x{aa,bb};
                    end
    
                end
                
           end
           
       end
        
    end
    
end



% Rounding relevant quantities
for aa=1:height(cases)
 
    cases.e(aa)=round(cases.e(aa),2);
    cases.dpdx{aa} = round(cases.dpdx{aa},0);
    cases.dpdx_exp{aa} = round(cases.dpdx_exp{aa},0);
    cases.mdot_f{aa} = round(cases.mdot_f{aa},3);
    cases.Q_f{aa} = round(cases.Q_f{aa},3);
    cases.U_f_x{aa} = round(cases.U_f_x{aa},3);
	cases.U_f_x_exp{aa} = round(cases.U_f_x_exp{aa},3);

end


save('wellboredata.mat', 'cases');



%% Add derived data

clear all;
close all;
clc;
load('wellboredata.mat', 'cases');

for aa=1:height(cases)
    % aa = 7;
    % aa = 35;
    % aa = 36;
     % aa = 47;
    % aa = 48;
    
    if strcmp(cases.flow_spec{aa},'U_f')
        cases.dpdx{aa} = zeros(length(cases.U_f_x{aa}),1);
        cases.Q_f{aa} = cases.U_f_x{aa}.*cases.A(aa);
        cases.mdot_f{aa} = cases.Q_f{aa}.*cases.rho_f(aa);
    end
    
    % Geometry coefficient beta
    if strcmp(cases.Geometry(aa),'Pipe')
        cases.beta(aa,1) = 3;
    else
        cases.beta(aa,1) = 2;
    end   
    
    
    
    % Loop rpm
    for bb = 1:length(cases.rpm{aa})   

        % Newtonian shear rate
        cases.SR_newt{aa,bb} = 24./cases.beta(aa).*cases.U_f_x{aa,bb}./cases.d_h(aa);

        % PL shear rate, handle n = 1 seperately due to formulation of PL shear
        % rate correction where for n_PL = 1 the denominator in the exponent
        % becomes zero
        if (cases.n_PL_f(aa)==1)
            cases.SR_PL{aa,bb} = cases.SR_newt{aa,bb};
        else
            cases.SR_PL{aa,bb} = cases.SR_newt{aa,bb} .* ((3.*cases.n_PL_f(aa)+1)./(4.*cases.n_PL_f(aa))) .^ (cases.n_PL_f(aa)./(cases.n_PL_f(aa)-1));
        end

        % PL viscosity
        cases.eta_PL{aa,bb} = cases.K_PL_f(aa).*cases.SR_PL{aa,bb}.^(cases.n_PL_f(aa)-1);

        % Generalized Reynolds number (Metzner-Reed 1955)
        % cases.Re_MR{aa,bb} = cases.rho_f(aa).*cases.U_f_x{aa,bb}.^(2-cases.n_PL_f(aa)).*cases.d_h(aa).^cases.n_PL_f(aa)./(cases.K_PL_f(aa).* ((3.*cases.n_PL_f(aa)+1)./(4.*cases.n_PL_f(aa))).^cases.n_PL_f(aa) .*8.^(cases.n_PL_f(aa)-1));

        % Generalized Reynolds number (Metzner-Reed 1955) with annular geometry correction of Kozicki et al. (1966) and Delplace and Leuliet (1995) 
        cases.Re_MR{aa,bb} = cases.rho_f(aa).*cases.U_f_x{aa,bb}.^(2-cases.n_PL_f(aa)).*cases.d_h(aa).^cases.n_PL_f(aa)./((24./cases.beta(aa)).^(cases.n_PL_f(aa)-1).*cases.K_PL_f(aa).*((1+cases.n_PL_f(aa).*cases.beta(aa))./(cases.n_PL_f(aa)+cases.n_PL_f(aa).*cases.beta(aa))).^cases.n_PL_f(aa));

        % Wall shear stress
        cases.tauw{aa,1} = cases.dpdx{aa,1}.*cases.d_h(aa)./4;

        % Dynamic pressure
        cases.p_dyn{aa,bb} = cases.rho_f(aa).*cases.U_f_x{aa,bb}.^2./2;

        % Friction factor
        cases.f{aa,bb} = cases.tauw{aa,1}./cases.p_dyn{aa,bb};

        % Friction velocity
        cases.ufricf{aa,bb} = sqrt(cases.tauw{aa,1}./cases.rho_f(aa));

        % Wall shear rate
        cases.SR_wall{aa,bb} = (cases.tauw{aa,1}./cases.K_PL_f(aa)).^(1./cases.n_PL_f(aa));

        % Wall viscosity
        cases.eta_wall{aa,bb} = cases.tauw{aa,1}./cases.SR_wall{aa,bb};

        % Reynolds number based on wall viscosity
        cases.Re_G{aa,bb} = cases.rho_f(aa).*cases.U_f_x{aa,bb}.*cases.d_h(aa)./cases.eta_wall{aa,bb};

    %     for bb=1:length(cases.ufricf{aa,bb})

            % u+
           %cases.u_plus{aa,bb}(bb,:)= cases.u_x{aa,bb}(bb,:)./cases.ufricf{aa,bb}(bb);

            % y+
            %cases.y_plus{aa,bb}(bb,:) = (cases.d_o(aa)/2 - cases.y{aa,bb}(bb,:)) .* cases.rho_f(aa).*cases.ufricf{aa,bb}(bb)./cases.eta_wall{aa,bb}(bb);

    %     end

        % Courant number 
        % cases.dt

        % ROP
        cases.ROP{aa,bb} = cases.A(aa).*cases.U_s_x{aa,bb} ./ (0.25.*pi.*cases.d_o(aa).^2) .* 3600;    
    
    end
    
end


% Rounding relevant quantities
for aa=1:height(cases)
 
    cases.e(aa)=round(cases.e(aa),2);
    cases.dpdx{aa} = round(cases.dpdx{aa},0);
    cases.dpdx_exp{aa} = round(cases.dpdx_exp{aa},0);
    cases.mdot_f{aa} = round(cases.mdot_f{aa},3);
    cases.Q_f{aa} = round(cases.Q_f{aa},3);
    cases.U_f_x{aa} = round(cases.U_f_x{aa},3);
	cases.U_f_x_exp{aa} = round(cases.U_f_x_exp{aa},3);

end


save('wellboredata.mat', 'cases');
