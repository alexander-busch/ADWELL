close all;
clear all;
clc;

addpath('utilities');
addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive';
filepath_metadata = [filepath '\metadata.xlsx'];
currentpath = pwd;


%% Import data and built table

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata);
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
cases = readtable(filepath_metadata,opts);


% Loop all table rows
for aa=1:height(cases)
    % Debugging 
    % aa = 4
    % aa = 5
    % aa = 31 
    % aa = height(cases)
    
    % Conversion from character strings to numerical arrays
    cases.dpdx{aa} = str2num(cases.dpdx{aa});
    cases.U_x{aa} = str2num(cases.U_x{aa});
    cases.rpm{aa} = str2num(cases.rpm{aa});       
    %cases.e{aa} = str2num(cases.e{aa}); 
    
    % compute x-sectional area
    cases.A{aa,1} = pi./4.*(cases.d_o(aa).^2-cases.d_i(aa).^2);
    
end

save('turbulencedata.mat', 'cases');

clear opts;


%% Set dpdx or mdot to zero depending on flow specification

clc;
load('turbulencedata.mat', 'cases');

% Loop all table rows
for aa=1:height(cases)
    % aa=height(cases)
    
    dpdx = cases.dpdx{aa};
    rpm = cases.rpm{aa};
    
    if strcmp(cases.Models{aa},'DNS')
        % DNS, dpdx is specified and mdot, U, Q are computed
        cases.mdot{aa,1} = [];
    else
        % CFD, either dpdx or mdto si specified depending on chosen
        % turbulence model
        if strcmp(cases.flow_spec{aa},'dpdx')
            cases.mdot{aa,1} = zeros(length(dpdx),1);
        else
            % Estimate mdot based on dpdx (Blasius)
            cases.U_x{aa,1} = (dpdx.^4.*cases.d_h(aa).^5./(1./16.*0.316.^4.*cases.K_PL_scaled(aa).*cases.rho_scaled(aa).^3)).^(1./7);
            cases.Q{aa,1} = cases.U_x{aa,1}.*cases.A{aa,1};
            cases.mdot{aa,1} = cases.Q{aa,1}.*cases.rho_scaled(aa);
            
            cases.dpdx{aa,1} = zeros(length(dpdx),1);

            % For some reason the following case has wrong mdots specified,
            % both in the filenames and in the actual CFD inputs...
            % aa = 9
            if strcmp(cases.ID(aa), '180526_3') 
                cases.mdot{aa,1} = [0.074946; 0.12652; 0.188; 0.27937; 0.32456; 0.34504; 0.39088; 0.41514; 0.52338; 0.6169];
            end

        end
    end 
end

save('turbulencedata.mat', 'cases');

clear aa data dpdx rpm;


%% Read simulation results

clc;
load('turbulencedata.mat', 'cases');

% Loop all table rows
for aa=1:height(cases)
    % aa=height(cases)
    % aa = 31;
    
    % Read data, check for DNS vs. CFD case
    if strcmp(cases.Turbulence(aa),'DNS')
        % DNS cases, read from metadata file
        cases.Q{aa,1} = cases.U_x{aa,1}.*cases.A{aa,1};
        cases.mdot{aa,1} = cases.Q{aa,1}.*cases.rho_scaled(aa);
    else
        % CFD cases
        
        dpdx = cases.dpdx{aa};
        rpm = cases.rpm{aa};
    
        % Determine sweeped quantity
        if length(dpdx)>length(rpm)
            sweep = dpdx;
        else
            sweep = rpm;
        end

        % Loop sweep
        for bb=1:length(sweep)
            % bb=2

            % Determine specification of dpdx or mdot and assemble respective filename
            if strcmp(cases.flow_spec{aa}, 'dpdx')
                if length(dpdx)>length(rpm)
                    filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_dpdx=' num2str(cases.dpdx{aa,1}(bb)) '_mass-flow-rate.txt'];
                else
                    filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_rpm=' num2str(cases.rpm{aa,1}(bb)) '_mass-flow-rate.txt'];
                end
            else    
                filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_mdot=' num2str(cases.mdot{aa,1}(bb)) '_dpdx.txt'];
            end

            % Check folder structure and assemble correct path to pressure loss/mass flow rate data
            if exist([filepath '\' cases.ID{aa} '\monitors'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\monitors'  filename];
            elseif exist([filepath '\' cases.ID{aa} '\' num2str(bb) '\monitors'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\' num2str(bb) '\monitors'  filename];
            else
                fullpath = [filepath '\' cases.ID{aa} '\' num2str(bb) '\solution\monitors'  filename];
            end

            % Check if file exists and read pressure loss/mass flow rate data
            if exist(fullpath,'file')
                if strcmp(cases.flow_spec{aa}, 'dpdx')                
                    cases.mdot{aa,1}(bb) = importMonitoredQuantities(fullpath);
                    
                    % Sanity check
                    if (bb>1) && (cases.mdot{aa,1}(bb)<cases.mdot{aa,1}(bb-1) && (length(dpdx)>length(rpm)) )
                        cases.mdot{aa,1}(1:bb-1)=NaN;
                    end
                    
                    if (strcmp(cases.ID{aa}, '180615_1') && bb==3)
                        cases.mdot{aa,1}(bb)=NaN;
                    end
                    
                    cases.Q{aa,1} = cases.mdot{aa,1}./cases.rho_scaled(aa);
                    cases.U_x{aa,1} = cases.Q{aa,1}./cases.A{aa,1};
                else 
                    cases.dpdx{aa,1}(bb) = importMonitoredQuantities(fullpath);

                    % Sanity check. For small Re_G, CFD simulations may have
                    % failed.
                    if ((cases.dpdx{aa,1}(bb)>0.1) ||  (cases.dpdx{aa,1}(bb)<0.0001))
                        cases.dpdx{aa,1}(bb)=NaN;
                        cases.mdot{aa,1}(bb)=NaN;
                    end
                end
            else
                cases.dpdx{aa,1}(bb)=0;
                cases.mdot{aa,1}(bb)=0;
                cases.Q{aa,1}=0;
                cases.U_x{aa,1}=0;
            end
        end
    end
end


% Loop all table rows and additionatilly read mdot_spec data
for aa=1:height(cases)
    % aa=height(cases)
    % aa = 9;
    
    fullpath = [filepath '\' cases.ID{aa} '\mdot_spec'];

    % Check if file exists and read pressure loss/mass flow rate data
    if exist(fullpath,'dir')
        
        % Estimate mdot based on dpdx (Blasius)
        U_x = (cases.dpdx{aa}.^4.*cases.d_h(aa).^5./(1./16.*0.316.^4.*cases.K_PL_scaled(aa).*cases.rho_scaled(aa).^3)).^(1./7);
        Q = U_x.*cases.A{aa,1};
        mdot = Q.*cases.rho_scaled(aa);
        
        % For some reason the following case has wrong mdots specified,
        % both in the filenames and in the actual CFD inputs...
        % aa = 9
        if strcmp(cases.ID(aa), '180526_3') 
            mdot = [0.074946; 0.12652; 0.188; 0.27937; 0.32456; 0.34504; 0.39088; 0.41514; 0.52338; 0.6169];
        end
        
        % Loop elements
        for bb=1:length(mdot)
            % bb=2

            filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_mdot=' num2str(mdot(bb)) '_dpdx.txt'];

            % Check folder structure and assemble correct path to pressure loss/mass flow rate data
            if exist([filepath '\' cases.ID{aa} '\mdot_spec\monitors'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\mdot_spec\monitors'  filename];
            elseif exist([filepath '\' cases.ID{aa} '\mdot_spec\' num2str(bb) '\monitors'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\mdot_spec\' num2str(bb) '\monitors'  filename];
            else
                fullpath = [filepath '\' cases.ID{aa} '\mdot_spec\' num2str(bb) '\solution\monitors'  filename];
            end
        
            % Check if file exists and read pressure loss/mass flow rate data
            if exist(fullpath,'file')
                dpdx(bb) = importMonitoredQuantities(fullpath);

                % Sanity check. For small Re_G, CFD simulations may have
                % failed.
                if ((cases.dpdx{aa,1}(bb)>0.1) ||  (cases.dpdx{aa,1}(bb)<0.0001))
                    dpdx(bb)=NaN;
                    mdot(bb)=NaN;
                end
            else
                cases.dpdx{aa,1}(bb)=0;
                cases.mdot{aa,1}(bb)=0;
                cases.Q{aa,1}=0;
                cases.U_x{aa,1}=0;
            end        
        end
        
        % Combine vectors
        %dpdx=dpdx';
        idx = mdot>min(cases.mdot{aa,1});
        mdot(idx)=[];
        dpdx(idx)=[];

        cases.dpdx{aa,1} = [dpdx; cases.dpdx{aa,1}];
        cases.mdot{aa,1} = [mdot; cases.mdot{aa,1}];
        
        cases.Q{aa,1} = cases.mdot{aa,1}./cases.rho_scaled(aa);
        cases.U_x{aa,1} = cases.Q{aa,1}./cases.A{aa,1};
        
    end
end
 

save('turbulencedata.mat', 'cases');

clear aa bb data dpdx rpm sweep;


%% Postprocess

clc;
load('turbulencedata.mat', 'cases');

% Loop all table rows
for aa=1:height(cases)
        
    % Geometry coefficient beta
    if strcmp(cases.Geometry(aa),'Pipe')
        cases.beta(aa) = 3;
    else
        cases.beta(aa) = 2;
    end   

    % Newtonian shear rate
    cases.SR_newt{aa,1} = 24./cases.beta(aa).*cases.U_x{aa,1}./cases.d_h(aa);
    
    % PL shear rate, handle n = 1 seperately due to formulation of PL shear
    % rate correction where for n_PL = 1 the denominator in the exponent
    % becomes zero
    if (cases.n_PL_scaled(aa)==1)
        cases.SR_PL{aa,1} = cases.SR_newt{aa,1};
    else
        cases.SR_PL{aa,1} = cases.SR_newt{aa,1} .* ((3.*cases.n_PL_scaled(aa)+1)./(4.*cases.n_PL_scaled(aa))) .^ (cases.n_PL_scaled(aa)./(cases.n_PL_scaled(aa)-1));
    end
    
    % PL viscosity
    cases.eta_PL{aa,1} = cases.K_PL_scaled(aa).*cases.SR_PL{aa,1}.^(cases.n_PL_scaled(aa)-1);
    
    % Generalized Reynolds number (Metzner-Reed 1955)
    % cases.Re_MR{aa,1} = cases.rho_scaled(aa).*cases.U_x{aa,1}.^(2-cases.n_PL_scaled(aa)).*cases.d_h(aa).^cases.n_PL_scaled(aa)./(cases.K_PL_scaled(aa).* ((3.*cases.n_PL_scaled(aa)+1)./(4.*cases.n_PL_scaled(aa))).^cases.n_PL_scaled(aa) .*8.^(cases.n_PL_scaled(aa)-1));
    
    % Generalized Reynolds number (Metzner-Reed 1955) with annular geometry correction of Kozicki et al. (1966) and Delplace and Leuliet (1995) 
    cases.Re_MR{aa,1} = cases.rho_scaled(aa).*cases.U_x{aa,1}.^(2-cases.n_PL_scaled(aa)).*cases.d_h(aa).^cases.n_PL_scaled(aa)./((24./cases.beta(aa)).^(cases.n_PL_scaled(aa)-1).*cases.K_PL_scaled(aa).*((1+cases.n_PL_scaled(aa).*cases.beta(aa))./(cases.n_PL_scaled(aa)+cases.n_PL_scaled(aa).*cases.beta(aa))).^cases.n_PL_scaled(aa));
    
     % Wall shear stress
    cases.tauw{aa,1} = cases.dpdx{aa,1}.*cases.d_h(aa)./4;

    % Dynamic pressure
    cases.p_dyn{aa,1} = cases.rho_scaled(aa).*cases.U_x{aa,1}.^2./2;
    
    % Friction factor
    cases.f{aa,1} = cases.tauw{aa,1}./cases.p_dyn{aa,1};
 
    % Friction velocity
    cases.ufricf{aa,1} = sqrt(cases.tauw{aa,1}./cases.rho_scaled(aa));

    % Wall shear rate
    cases.SR_wall{aa,1} = (cases.tauw{aa,1}./cases.K_PL_scaled(aa)).^(1./cases.n_PL_scaled(aa));

    % Wall viscosity
    cases.eta_wall{aa,1} = cases.tauw{aa,1}./cases.SR_wall{aa,1};
    
    % Reynolds number based on viscosity at the wall
    cases.Re_G{aa,1} = cases.rho_scaled(aa).*cases.U_x{aa,1}.*cases.d_h(aa)./cases.eta_wall{aa,1};

    
        % Sanity check
    if (~strcmp(cases.Turbulence{aa},'DNS') &&length(cases.rpm{aa})==1 && ~issorted(cases.f{aa,1}))
        for bb=2:length(cases.dpdx{aa})
           if cases.f{aa}(bb)>cases.f{aa}(bb-1)
               cases.f{aa}(bb-1)=NaN;
               cases.Re_MR{aa}(bb-1)=NaN;
               cases.Re_G{aa}(bb-1)=NaN;
               cases.mdot{aa}(bb-1)=NaN;
               cases.U_x{aa}(bb-1)=NaN;
               cases.Q{aa}(bb-1)=NaN;
           end
               
               
        end
        
    end

end

save('turbulencedata.mat', 'cases');

clear aa bb a b data dpdx rpm sweep;



%% Read CFD velocity profile data and compute u+ and y+

clc;
load('turbulencedata.mat', 'cases');

% figure 
% hold on

% Loop all table rows
for aa=1:height(cases)
    % aa = 20
    dpdx = cases.dpdx{aa};
    rpm = cases.rpm{aa};

    % Determine sweeped quantity
    if length(dpdx)>length(rpm)
        sweep = dpdx;
    else
        sweep = rpm;
    end

    % Create table variables for y-coordinate and x-velocity as well as y+ and u+ this will
    % hold arrays corresponding to dpdx
    cases.y{aa,1}=zeros(length(dpdx),100);
    cases.u_x{aa,1}=zeros(length(dpdx),100);
	cases.y_plus{aa,1}=zeros(length(dpdx),100);
    cases.u_plus{aa,1}=zeros(length(dpdx),100);

    % Loop sweep and get data
    for bb=1:length(sweep)
        % bb = 6
        if strcmp(cases.ID(aa),'180521_1') %     bb~=3
            % Skip as this case has buggy fluent export files 
        else

            % Determine specification of dpdx or mdot and assemble respective filename
            if strcmp(cases.flow_spec{aa}, 'dpdx')
                if length(dpdx)>length(rpm)
                    filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_dpdx=' num2str(cases.dpdx{aa,1}(bb)) '.txt'];
                else
                    filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_rpm=' num2str(cases.rpm{aa,1}(bb)) '.txt'];
                end
            else    
                filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_mdot=' num2str(cases.mdot{aa,1}(bb)) '.txt'];
            end

            % Check folder structure and assemble correct path to pressure loss/mass flow rate data
            if exist([filepath '\' cases.ID{aa} '\exports'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\exports'  filename];
            elseif exist([filepath '\' cases.ID{aa} '\' num2str(bb) '\exports'],'dir')
                fullpath = [filepath '\' cases.ID{aa} '\' num2str(bb) '\exports'  filename];
            else
                fullpath = [filepath '\' cases.ID{aa} '\' num2str(bb) '\solution\exports'  filename];
            end

            % Check if file exists and read velocity profiles
            if exist(fullpath,'file')
                % [cases.y{aa,1}(bb,:),cases.u_x{aa,1}(bb,:)] = importExportedQuantities(fullpath);
                [y, u_x] = importExportedQuantities(fullpath, cases.ID{aa});
                
                % Remove NaNs
                y = (y(~isnan(y)));
                u_x = (u_x(~isnan(u_x)));

                % Flip matrices in case of wrong ordering
                if y(1)>y(2)
                   y = fliplr(y);
                   u_x = fliplr(u_x);
                end
                
                % Add missing data point for u_x at y=0 in case of pipe
                if (strcmp(cases.Geometry{aa},'Pipe') && ~all(u_x))
  
                    % Interpolate u_x at y = 0 and add to data
                    u_x = [u_x(sign(y) < 0) interp1(y(~isnan(y)),u_x(~isnan(u_x)),0,'spline') u_x(sign(y) > 0)];
                    y = [y(sign(y) < 0) 0 y(sign(y) > 0)];
                
                end

                % Write to table
                cases.y{aa,1}(bb,:) = NaN;
                cases.u_x{aa,1}(bb,:) = NaN;
                cases.y{aa,1}(bb,1:length(y)) = y;
                cases.u_x{aa,1}(bb,1:length(u_x)) = u_x;
                
                clear y u_x;                    


%                 Sanity check. For small Re_G, CFD simulations may have
%                 failed.
%                 if ((cases.dpdx{aa,1}(bb)>0.1) ||  (cases.dpdx{aa,1}(bb)<0.001))
%                     cases.y{aa,1}(bb)=0;
%                     cases.u_x{aa,1}(bb)=0;
%                 end

            else
                cases.y{aa,1}(bb,:)=0;
                cases.u_x{aa,1}(bb,:)=0;
            end
            
        % Uncomment ot plot all u_x y data on on single figure
%         scatter(cases.y{aa,1}(bb,:),cases.u_x{aa,1}(bb,:))
       
        end
    end
  
    % u+ and y+
    for bb=1:length(cases.ufricf{aa,1})
        
        % Determine u+ and y+
        cases.u_plus{aa,1}(bb,:)= cases.u_x{aa,1}(bb,:)./cases.ufricf{aa,1}(bb);
        %cases.y_plus{aa,1}(bb,:) = (cases.d_o(aa)/2 - cases.y{aa,1}(bb,:)) .* cases.rho_scaled(aa).*cases.ufricf{aa,1}(bb)./cases.eta_wall{aa,1}(bb);
        cases.y_plus{aa,1}(bb,:) = (cases.d_o(aa)/2 + cases.y{aa,1}(bb,:)) .* cases.rho_scaled(aa).*cases.ufricf{aa,1}(bb)./cases.eta_wall{aa,1}(bb);
        
        Rplus = cases.rho_scaled(aa).*cases.ufricf{aa,1}(bb)./cases.eta_wall{aa,1}(bb)/ 2;
            
        % Cut off at max(uplus)
        idx = find(cases.u_plus{aa,1}(bb,:) == max(cases.u_plus{aa,1}(bb,:)));
        cases.u_plus{aa,1}(bb,idx:end) = NaN;
        cases.y_plus{aa,1}(bb,idx:end) = NaN;           

        % Uncomment ot plot all uplus yplus data on on single figure
        % scatter(cases.y_plus{aa,1}(bb,:),cases.u_plus{aa,1}(bb,:))
    end
    
end

% aa=18
% bb=6
% scatter(cases.y{aa,1}(bb,:),cases.u_x{aa,1}(bb,:))
% scatter(cases.y_plus{aa,1}(bb,:),cases.u_plus{aa,1}(bb,:))


save('turbulencedata.mat', 'cases');

clear aa bb a b data dpdx rpm sweep;


%% Extract data from DNS velocity profile and u+ y+ figure

clc;
load('turbulencedata.mat', 'cases');

% h = open('C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive\DNS velocity profiles\Comp_fz0144_yUz.fig');
h = open('C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive\DNS velocity profiles\Comp_fz0144_yUz_corrected.fig');

axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axe
dataObjs = dataObjs{2};
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

% DNS velocity profiles are normalized with bulk velocity (velocity) and
% pipe radius (y-coordinate), hence DNS profiles are multiplied with the
% respective quantities

% PL
aa=15;
bb=6;
cases.y{aa,1}(bb,1:length(cell2mat(xdata(2)))) = cell2mat(xdata(2))*(cases.d_o(aa)/2);
cases.u_x{aa,1}(bb,1:length(cell2mat(ydata(2)))) = cell2mat(ydata(2))*cases.U_x{aa}(bb);

% Cross
aa=16;
bb=6;
cases.y{aa,1}(bb,1:length(cell2mat(xdata(1)))) = cell2mat(xdata(1))*(cases.d_o(aa)/2);
cases.u_x{aa,1}(bb,1:length(cell2mat(ydata(1)))) = cell2mat(ydata(1))*cases.U_x{aa}(bb);

% Newtonian
aa=17;
bb=6;
cases.y{aa,1}(bb,1:length(cell2mat(xdata(3)))) = cell2mat(xdata(3))*(cases.d_o(aa)/2);
cases.u_x{aa,1}(bb,1:length(cell2mat(ydata(3)))) = cell2mat(ydata(3))*cases.U_x{aa}(bb);
%  index (3) is the corrected Newtonian case, index (4) is the
%  old/incorrect Newtonain case
% test = figure();
% hold on;
% grid on;
% for ii = 1:3
%     % ii = 4
%     plot(cell2mat(xdata(ii)),cell2mat(ydata(ii)));
% end

close(h);

h = open('C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive\DNS velocity profiles\Comp_fz0144_ypUp.fig');
set(gca, 'xlim', [0 300],...
    'XScale', 'log',...  
    'YScale', 'lin');
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axe
dataObjs = dataObjs{2};
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

aa=15;
bb=6;
cases.y_plus{aa,1}(bb,1:length(cell2mat(xdata(2)))) = cell2mat(xdata(2));
cases.u_plus{aa,1}(bb,1:length(cell2mat(ydata(2)))) = cell2mat(ydata(2));

aa=16;
bb=6;
cases.y_plus{aa,1}(bb,1:length(cell2mat(xdata(1)))) = cell2mat(xdata(1));
cases.u_plus{aa,1}(bb,1:length(cell2mat(ydata(1)))) = cell2mat(ydata(1));

aa=17;
bb=6;
cases.y_plus{aa,1}(bb,1:length(cell2mat(xdata(3)))) = cell2mat(xdata(3));
cases.u_plus{aa,1}(bb,1:length(cell2mat(ydata(3)))) = cell2mat(ydata(3));

% viscsublay = [ cell2mat(xdata(5))', cell2mat(ydata(5))'];
% loglaw = [ cell2mat(xdata(4))', cell2mat(ydata(4))'];

close(h);

save('turbulencedata.mat', 'cases');

clear aa bb h axesObjs dataObjs objTypes xdata ydata;





