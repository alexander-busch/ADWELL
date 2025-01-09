% Create hpc array jobs journal files for turbulent pipe and annular flow
% Modify xlsread(filepath_metadata,'CFD','DC135:DH144'); to read correct cells
% Modify absolute paths of filepath_metadata & filepath_journal

clear all;
close all;
clc;

filepath_metadata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive\metadata.xlsx';

% Define case to simulate
case2sim = '181030_1';

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata);
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
cases = readtable(filepath_metadata,opts);

% Conversion from character strings to numerical arrays
for aa=1:height(cases)
    cases.rpm{aa} = str2num(cases.rpm{aa});
    cases.dpdx{aa} = str2num(cases.dpdx{aa});
    cases.U_x{aa} = str2num(cases.U_x{aa});
end


% Select 
case2sim = cases(strcmp(cases.ID,case2sim),:);

% Parameters
if strcmp(case2sim.Geometry,'Pipe')
    geo = 'p'; 
    dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04]; % dpdx for pipe DNS simulations
else
    geo = 'a';

    if length(case2sim.dpdx{:})==1
        sweep = case2sim.rpm{:};
        dpdx = case2sim.dpdx{:};
    else
        dpdx = case2sim.dpdx{:};
        sweep = case2sim.dpdx{:};
    end       
        
end

if contains(case2sim.Turbulence,'SST')
    turbmod = 'omega'; 
else
    turbmod = 'epsilon';
end

fluidlist = unique(cases.Fluid);
fluidnames = {'cross',...
    'h2o',...
    'newt',...
    'pl'};
for ii=1:length(fluidlist)
    if strcmp(case2sim.Fluid{:},fluidlist{ii})
        fluid = fluidnames{ii};
    end
end

% Estimate mdot based on dpdx (Blasius)
vel = (dpdx.^4.*case2sim.d_h.^5./(1./16.*0.316.^4.*case2sim.K_PL_scaled.*case2sim.rho_scaled.^3)).^(1./7);
A = pi./4.*(case2sim.d_o.^2-case2sim.d_i.^2);
Q = vel.*A;
mdot = Q.*case2sim.rho_scaled;

% Path and filename of journalfile
filepath_journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\hpc\journals_hpc';

% Purge
if exist(filepath_journal,'dir')
    rmdir(filepath_journal,'s');
end
mkdir(filepath_journal);

% Loop all dpdx
for ii = 1:length(sweep)

    fileID = fopen([filepath_journal '\run' num2str(ii) '.jou'],'w');

    % Write journal file
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Pressure gradient or rpm sweep of turbulent pipe/annular flow\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; 02.10.2018, alexander.busch@alumni.ntnu.no\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; Reads generic case and sweeps pressure gradient or rpm\r\n');
    fprintf(fileID, '; Attention: Change absolute path in rename command to suit working directory \r\n');
    fprintf(fileID, ['; Copy journal file in Fluent working directory and paste "file/read-journal run' num2str(ii) '.jou" to run simulations for defined pressure gradients/mass flow rates\r\n']);
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Read case and change color scheme\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
%     if strcmp(geo,'p')
%         fprintf(fileID, '/file/read-case-data pipe.cas OK\r\n');
%     else
%         fprintf(fileID, '/file/read-case-data annulus.cas OK\r\n'); 
%     end
    fprintf(fileID, '/file/read-case case.cas\r\n'); 
    fprintf(fileID, '/display/set/colors/color-scheme classic\r\n');
    fprintf(fileID, '/display/set/colors/background "white"\r\n');
    fprintf(fileID, '/display/set/colors/foreground "black"\r\n');
    fprintf(fileID, '/display/re-render\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Purge working directory\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '!if exist solution rmdir /s /q solution\r\n');
    fprintf(fileID, '!mkdir solution\r\n');
    fprintf(fileID, '!mkdir "solution/reports"\r\n');
    fprintf(fileID, '!mkdir "solution/exports"\r\n');
    fprintf(fileID, '!mkdir "solution/monitors"\r\n');
    fprintf(fileID, '!mkdir "solution/autosave"\r\n');
    fprintf(fileID, '!mkdir "solution/cfd-post"\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Parameters\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    if strcmp(fluid,'newt')
        fprintf(fileID, ['(rpsetvar ''newt/eta/scaled ' num2str(case2sim.K_PL_scaled(:)) ')\r\n']);
        fprintf(fileID, '/define/materials/change-create newt-scaled newt-scaled y constant (%%rpgetvar ''newt/rho/scaled) n n y constant (%%rpgetvar ''newt/eta/scaled) n n n\r\n');
    elseif strcmp(fluid,'pl')
        
        
    elseif strcmp(fluid,'cross')
        fprintf(fileID, '(if (not (rp-var-object ''cross/k/scaled))		(rp-var-define ''cross/k/scaled 9.58246 ''real #f))\r\n');
        fprintf(fileID, '(if (not (rp-var-object ''cross/n/scaled))		(rp-var-define ''cross/n/scaled 0.60888 ''real #f))\r\n');
        fprintf(fileID, '(if (not (rp-var-object ''cross/mu_0/scaled))	(rp-var-define ''cross/mu_0/scaled 9.4807e-04 ''real #f))\r\n');
        fprintf(fileID, '(if (not (rp-var-object ''cross/mu_inf/scaled))	(rp-var-define ''cross/mu_inf/scaled 9.4807e-06 ''real #f))\r\n');
    else
        
    end

    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, ['; Solve - Scaled ' fluid ' material function\r\n']);
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');

    if strcmp(fluid,'cross')
        fprintf(fileID, '/define/user-defined/compiled-functions load "libudf"\r\n');
    end
    fprintf(fileID, ['/define/boundary-conditions/fluid fluid y ' fluid '-scaled , , , , , , , , , , , , , , , , , , , ,\r\n']);
    fprintf(fileID, '\r\n');
    
    if length(case2sim.dpdx{:})>1

        if strcmp(turbmod,'epsilon')
            fprintf(fileID, ['/define/periodic-conditions/massflow-rate-specification ' num2str(m_dot(ii)) ' -' num2str(dpdx(ii)) ' 0.5 5 1 0 0\r\n']);
            fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 3 y y "solution/monitors/' geo '_' fluid '-scaled_mdot=' num2str(m_dot(ii)) '_mass-flow-rate.txt" 1\r\n']);
            fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 4 y y "solution/monitors/' geo '_' fluid '-scaled_mdot=' num2str(m_dot(ii)) '_bulk-velocity.txt" 1\r\n']);
            fprintf(fileID, '/solve/monitors/statistic/monitors y\r\n'); 
            fprintf(fileID, '/solve/monitors/statistic/print y\r\n'); 
            fprintf(fileID, '/solve/monitors/statistic/plot y\r\n'); 
            fprintf(fileID, '/solve/monitors/statistic/write y\r\n'); 
            fprintf(fileID, '/solve/monitors/statistic/window 5\r\n');
        else
            fprintf(fileID, ['/define/periodic-conditions/pressure-gradient-specification -' num2str(dpdx(ii)) ' 1 0 0\r\n']);
            fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 2 y y "solution/monitors/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '_mass-flow-rate.txt" 1\r\n']);
            fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 3 y y "solution/monitors/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '_bulk-velocity.txt" 1\r\n']);
        end
        
            fprintf(fileID, ['/solve/initialize/set-defaults/x-velocity ' num2str(vel(ii)) '\r\n']);
        
    else
        % Sweep rpm
        
        fprintf(fileID, ['/define/periodic-conditions/pressure-gradient-specification -' num2str(dpdx) ' 1 0 0\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 2 y y "solution/monitors/' geo '_' fluid '-scaled_rpm=' num2str(sweep(ii)) '_mass-flow-rate.txt" 1\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 3 y y "solution/monitors/' geo '_' fluid '-scaled_rpm=' num2str(sweep(ii)) '_bulk-velocity.txt" 1\r\n']);
        fprintf(fileID, ['/define/boundary-conditions/set/wall wall_inner () motion-bc y motion-bc-moving rotating y omega n ' num2str(2*pi*sweep(ii)/60) ' ai 1 aj 0 ak 0 q\r\n']);
        
            fprintf(fileID, ['/solve/initialize/set-defaults/x-velocity ' num2str(vel) '\r\n']);
    end
    

    fprintf(fileID, '/solve/initialize/set-defaults/k 1\r\n'); % num2str(k(ii)) '\r\n']);
    fprintf(fileID, ['/solve/initialize/set-defaults/' turbmod ' 1 \r\n']); % ' num2str(epsilon(ii)) '\r\n']);
    fprintf(fileID, '/solve/initialize/initialize-flow ok ,\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme pressure 10\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme mom 0\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme k 0\r\n');
%     fprintf(fileID, ['/solve/set/discretization-scheme ' turbmod ' 0\r\n']);
%     fprintf(fileID, '/solve/iterate 4000\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme pressure 12\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme mom 1\r\n');
%     fprintf(fileID, '/solve/set/discretization-scheme k 1\r\n');
%     fprintf(fileID, ['/solve/set/discretization-scheme ' turbmod ' 1\r\n']);
    fprintf(fileID, ['/file/write-case-data "case_init.cas"\r\n']);    
    fprintf(fileID, '/solve/iterate 20000\r\n');

    if length(case2sim.dpdx{:})>1
        
        if strcmp(turbmod,'epsilon')
            fprintf(fileID, ['/file/export/ascii "solution/exports/' geo '_' fluid '-scaled_m-dot=' num2str(m_dot(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
            fprintf(fileID, ['/report/summary y "solution/reports/' geo '_' fluid '-scaled_m-dot=' num2str(m_dot(ii)) '.txt"\r\n']);
            fprintf(fileID, ['/file/write-case-data "solution/' geo '_' fluid '-scaled_m-dot=' num2str(m_dot(ii)) '.cas"\r\n']);
            fprintf(fileID, '!move /y statistics-per-pr-grad.out solution\\monitors\\statistics-per-pr-grad.out\r\n');
            fprintf(fileID, ['!rename ' strrep(filepath_journal,'\','\\') '\\solution\\monitors\\statistics-per-pr-grad.out "' geo '_' fluid '-scaled_mdot=' num2str(m_dot(ii)) '_dpdx.txt"\r\n']);
        else
            fprintf(fileID, ['/file/export/ascii "solution/exports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate specific-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
            fprintf(fileID, ['/report/summary y "solution/reports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt"\r\n']);
            fprintf(fileID, ['/file/write-case-data "solution/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.cas"\r\n']);
        end
    
    else
            fprintf(fileID, ['/file/export/ascii "solution/exports/' geo '_' fluid '-scaled_rpm=' num2str(sweep(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate specific-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
            fprintf(fileID, ['/report/summary y "solution/reports/' geo '_' fluid '-scaled_rpm=' num2str(sweep(ii)) '.txt"\r\n']);
            fprintf(fileID, ['/file/write-case-data "solution/' geo '_' fluid '-scaled_rpm=' num2str(sweep(ii)) '.cas"\r\n']);
    end
    
    fprintf(fileID, '\r\n');
    
    fclose(fileID);
end

