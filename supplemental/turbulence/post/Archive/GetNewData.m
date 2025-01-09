close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive';
filepath_metadata = [filepath '\metadata.xlsx'];
currentpath = pwd;


%% Import data and built table

% Read existing data
load('turbulencedata.mat');

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata);
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
newcases = readtable(filepath_metadata,opts);

if height(cases)<height(newcases)
    newcases(1:height(cases),:) = [];
end
% cd(filepath); 
% cd(currentpath); 

% Loop all new table rows
for aa=1:height(newcases)
    % Debugging 
    % aa = 4
    % aa = 5
    % aa = 40
    
    % Conversion from character strings to numerical arrays
    newcases.dpdx{aa} = str2num(newcases.dpdx{aa});
    newcases.U_x{aa} = str2num(newcases.U_x{aa});
    newcases.rpm{aa} = str2num(newcases.rpm{aa});   
    
    % compute x-sectional area
    newcases.A{aa,1} = pi./4.*(newcases.d_o(aa).^2-newcases.d_i(aa).^2);
    
    % Get dpdx write estimated U and mdot to table
    
%     if (strcmp(newcases.Models{aa},'DNS') && strcmp(newcases.Geometry{aa},'Annulus'))
%         newcases.dpdx{aa,1} = newcases.dpdx_DNS(aa,1);
%     else
%         if strcmp(newcases.Geometry{aa},'Pipe')
%             % Predefined pressure gradient [Pa] used in all pipe DNS and k-omega CFD
%             % simulations. In case of the k-epsilon simulations, a mass flow
%             % rate estimated based on dpdx (Blasius 1912) was used as input instead.
%             dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04]; % dpdx for pipe DNS simulations
%         else
%             if contains(newcases.ID{aa},'180605')
%                 % Old annular cases were run with same dpdx as pipe cases
%                 dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04]; % dpdx for pipe DNS simulations
%             else               
%                 % Predefined pressure gradient [Pa] used in all Newtonian k-omega CFD
%                 % simulations (based on the annular Newtonian DNS and extended to lower Re).
%                 dpdx = [0.00207; 0.00414; 0.00828; 0.01656]; % dpdx for annular DNS simulations 
%                 dpdx = [0.00075; 0.00075; 0.001; dpdx ]; % Enlarged for lower Re
%             end
%         end
%         
%         
%         newcases.dpdx{aa,1} = dpdx;
%         
%     end

    dpdx = cases.dpdx{aa};
    newcases.mdot{aa,1} = zeros(length(dpdx),1);
    
    if strcmp(newcases.Models{aa},'DNS') 
        newcases.mdot{aa,1} = [];
    else
    
        % Estimate mdot based on dpdx (Blasius)
        newcases.U_x{aa,1} = (dpdx.^4.*newcases.d_h(aa).^5./(1./16.*0.316.^4.*newcases.K_PL_scaled(aa).*newcases.rho_scaled(aa).^3)).^(1./7);
        newcases.Q{aa,1} = newcases.U_x{aa,1}.*newcases.A{aa,1};
        newcases.mdot{aa,1} = newcases.Q{aa,1}.*newcases.rho_scaled(aa);

        % Overwrite with zeros to represent specification
        if strcmp(newcases.flow_spec{aa},'dpdx')
            newcases.mdot{aa,1} = zeros(length(dpdx),1);
        else
            newcases.dpdx{aa,1} = zeros(length(dpdx),1);

            % For some reason the following case has wrong mdots specified,
            % both in the filenames and in the actual CFD inputs...
            if strcmp(newcases.ID(aa), '180526_3') 
                newcases.mdot{aa,1} = [0.074946; 0.12652; 0.188; 0.27937; 0.32456; 0.34504; 0.39088; 0.41514; 0.52338; 0.6169];
            end
        end
    end
    
    
    % Create table variables for y-coordinate and x-velocity, this will
    % hold arrays corresponding to dpdx
    newcases.y{aa,1}=zeros(length(dpdx),100);
    newcases.u_x{aa,1}=zeros(length(dpdx),100);
    
    
    % Read data, check for DNS vs. CFD case
    if strcmp(newcases.Turbulence(aa),'DNS') % DNS newcases, read from metadata file
        if strcmp(newcases.Geometry{aa},'Pipe')
            data = xlsread(filepath_metadata,newcases.ID{aa},'A2:B11');
            newcases.dpdx{aa,1}=data(:,1);
            newcases.mdot{aa,1}=data(:,2);
        else 
            newcases.U_x{aa,1} = newcases.U_x_mean_DNS(aa);
            newcases.Q{aa,1} = newcases.U_x{aa,1}.*newcases.A{aa,1};
            newcases.mdot{aa,1} = newcases.Q{aa,1}.*newcases.rho_scaled(aa);
        end
    else % CFD newcases
        
        % Loop all dpdx
        for bb=1:length(dpdx)
            
            % Pressure loss/mass flow rate data ---------------------------

            % Determine specification of dpdx or mdot and assemble respective filename
            if strcmp(newcases.flow_spec{aa}, 'dpdx')
                filename = ['\' newcases.Geometry{aa}(1) '_' newcases.Fluid{aa} '-scaled_dpdx=' num2str(newcases.dpdx{aa,1}(bb)) '_mass-flow-rate.txt'];
            else    
                filename = ['\' newcases.Geometry{aa}(1) '_' newcases.Fluid{aa} '-scaled_mdot=' num2str(newcases.mdot{aa,1}(bb)) '_dpdx.txt'];
            end

            % Change to directory 
            % cd([filepath '\' newcases.ID{aa}]);

            % Check folder structure and assemble correct path to pressure loss/mass flow rate data
            if exist([filepath '\' newcases.ID{aa} '\monitors'],'dir')
                fullpath = [filepath '\' newcases.ID{aa} '\monitors'  filename];
            elseif exist([filepath '\' newcases.ID{aa} '\' num2str(bb) '\monitors'],'dir')
                fullpath = [filepath '\' newcases.ID{aa} '\' num2str(bb) '\monitors'  filename];
            else
                fullpath = [filepath '\' newcases.ID{aa} '\' num2str(bb) '\solution\monitors'  filename];
            end

            % Check if file exists and read pressure loss/mass flow rate data
            if exist(fullpath,'file')
                if strcmp(newcases.flow_spec{aa}, 'dpdx')                
                    newcases.mdot{aa,1}(bb) = importMonitoredQuantities(fullpath);
                else 
                    newcases.dpdx{aa,1}(bb) = importMonitoredQuantities(fullpath);

                    % Sanity check. For small Re_G, CFD simulations may have
                    % failed.
%                     if ((newnewcases.dpdx{aa,1}(bb)>0.1) ||  (newnewcases.dpdx{aa,1}(bb)<0.001))
%                         newnewcases.dpdx{aa,1}(bb)=0;
%                         newnewcases.mdot{aa,1}(bb)=0;
%                     end
                end
            else
                newcases.dpdx{aa,1}(bb)=0;
                newcases.mdot{aa,1}(bb)=0;
            end
            
            
            % Velocity profile data ---------------------------------------
            
            if strcmp(newcases.ID(aa),'180521_1') %     bb~=3
                % Skip as this case has buggy fluent export files 
            else
            
                % Determine specification of dpdx or mdot and assemble respective filename
                if strcmp(newcases.flow_spec{aa}, 'dpdx')
                    filename = ['\' newcases.Geometry{aa}(1) '_' newcases.Fluid{aa} '-scaled_dpdx=' num2str(newcases.dpdx{aa,1}(bb)) '.txt'];
                else    
                    filename = ['\' newcases.Geometry{aa}(1) '_' newcases.Fluid{aa} '-scaled_mdot=' num2str(newcases.mdot{aa,1}(bb)) '.txt'];
                end

                % Check folder structure and assemble correct path to pressure loss/mass flow rate data
                if exist([filepath '\' newcases.ID{aa} '\exports'],'dir')
                    fullpath = [filepath '\' newcases.ID{aa} '\exports'  filename];
                elseif exist([filepath '\' newcases.ID{aa} '\' num2str(bb) '\exports'],'dir')
                    fullpath = [filepath '\' newcases.ID{aa} '\' num2str(bb) '\exports'  filename];
                else
                    fullpath = [filepath '\' newcases.ID{aa} '\' num2str(bb) '\solution\exports'  filename];
                end

                % Check if file exists and read velocity profiles
                if exist(fullpath,'file')
                    % [cases.y{aa,1}(bb,:),cases.u_x{aa,1}(bb,:)] = importExportedQuantities(fullpath);
                    [a, b] = importExportedQuantities(fullpath);
                    newcases.y{aa,1}(bb,1:length(a)) = a;
                    newcases.u_x{aa,1}(bb,1:length(b)) = b;
                    clear a b;

                    % Sanity check. For small Re_G, CFD simulations may have
                    % failed.
    %                 if ((newcases.dpdx{aa,1}(bb)>0.1) ||  (newcases.dpdx{aa,1}(bb)<0.001))
    %                     newcases.dpdx{aa,1}(bb)=0;
    %                     newcases.mdot{aa,1}(bb)=0;
    %                 end

                else
                    newcases.y{aa,1}(bb,:)=0;
                    newcases.u_x{aa,1}(bb,:)=0;
                end
            end
        end
    end
end

% aa=5
% bb=2
% plot(cases.y{aa,1}(bb,:),newcases.u_x{aa,1}(bb))


%% Postprocess

% Create table variables for y-coordinate and x-velocity, this will
% hold arrays corresponding to dpdx
newcases.u_plus{aa,1}=zeros(length(dpdx),100);
newcases.y_plus{aa,1}=zeros(length(dpdx),100);

% Loop all table rows
for aa=1:height(newcases)
    % Recalculate volumetric flow rate and bulk velocity
    newcases.Q{aa,1} = newcases.mdot{aa,1}./newcases.rho_scaled(aa);
    newcases.U_x{aa,1} = newcases.Q{aa,1}./newcases.A{aa,1};
        
    % Geometry coefficient beta
    if strcmp(newcases.Geometry(aa),'Pipe')
        newcases.beta(aa) = 3;
    else
        newcases.beta(aa) = 2;
    end   

    % Newtonian shear rate
    newcases.SR_newt{aa,1} = 24./newcases.beta(aa).*newcases.U_x{aa,1}./newcases.d_h(aa);
    
    % PL shear rate, handle n = 1 seperately due to formulation of PL shear
    % rate correction where for n_PL = 1 the denominator in the exponent
    % becomes zero
    if (newcases.n_PL_scaled(aa)==1)
        newcases.SR_PL{aa,1} = newcases.SR_newt{aa,1};
    else
        newcases.SR_PL{aa,1} = newcases.SR_newt{aa,1} .* ((3.*newcases.n_PL_scaled(aa)+1)./(4.*newcases.n_PL_scaled(aa))) .^ (newcases.n_PL_scaled(aa)./(newcases.n_PL_scaled(aa)-1));
    end
    
    % PL viscosity
    newcases.eta_PL{aa,1} = newcases.K_PL_scaled(aa).*newcases.SR_PL{aa,1}.^(newcases.n_PL_scaled(aa)-1);
    
    % Reynolds number based on apparent viscosity
    % newcases.Re_G{aa,1} = newcases.rho_scaled(aa).*newcases.U_x{aa,1}.*newcases.d_h(aa)./newcases.eta_PL{aa,1};
    
    % Generalized Reynolds number (Metzner-Reed 1955)
    % newcases.Re_G{aa,1} = newcases.rho_scaled(aa).*newcases.U_x{aa,1}.^(2-newcases.n_PL_scaled(aa)).*newcases.d_h(aa).^newcases.n_PL_scaled(aa)./(newcases.K_PL_scaled(aa).* ((3.*newcases.n_PL_scaled(aa)+1)./(4.*newcases.n_PL_scaled(aa))).^newcases.n_PL_scaled(aa) .*8.^(newcases.n_PL_scaled(aa)-1));
    
    % Generalized Reynolds number (Metzner-Reed 1955) with annular geometry correction of Kozicki et al. (1966) and Delplace and Leuliet (1995) 
    newcases.Re_G{aa,1} = newcases.rho_scaled(aa).*newcases.U_x{aa,1}.^(2-newcases.n_PL_scaled(aa)).*newcases.d_h(aa).^newcases.n_PL_scaled(aa)./((24./newcases.beta(aa)).^(newcases.n_PL_scaled(aa)-1).*newcases.K_PL_scaled(aa).*((1+newcases.n_PL_scaled(aa).*newcases.beta(aa))./(newcases.n_PL_scaled(aa)+newcases.n_PL_scaled(aa).*newcases.beta(aa))).^newcases.n_PL_scaled(aa));
    
     % Wall shear stress
    newcases.tauw{aa,1} = newcases.dpdx{aa,1}.*newcases.d_h(aa)./4;

    % Dynamic pressure
    newcases.p_dyn{aa,1} = newcases.rho_scaled(aa).*newcases.U_x{aa,1}.^2./2;
    
    % Friction factor
    newcases.f{aa,1} = newcases.tauw{aa,1}./newcases.p_dyn{aa,1};

    % Friction velocity
    newcases.ufricf{aa,1} = sqrt(newcases.tauw{aa,1}./newcases.rho_scaled(aa));

    % Wall shear rate
    newcases.SR_wall{aa,1} = (newcases.tauw{aa,1}./newcases.K_PL_scaled(aa)).^(1./newcases.n_PL_scaled(aa));

    % Wall viscosity
    newcases.eta_wall{aa,1} = newcases.tauw{aa,1}./newcases.SR_wall{aa,1};


    for bb=1:length(newcases.ufricf{aa,1})

        % u+
        newcases.u_plus{aa,1}(bb,:)= newcases.u_x{aa,1}(bb,:)./newcases.ufricf{aa,1}(bb);

        % y+
        newcases.y_plus{aa,1}(bb,:) = (newcases.d_o(aa)/2 - newcases.y{aa,1}(bb,:)) .* newcases.rho_scaled(aa).*newcases.ufricf{aa,1}(bb)./newcases.eta_wall{aa,1}(bb);

    end

end


%% Assemble final table 

cases = [cases; newcases];

save('turbulencedata.mat', 'cases');