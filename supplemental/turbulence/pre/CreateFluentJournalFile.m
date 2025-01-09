% If you plan to read the file with Microsoft® Notepad, use '\r\n' instead of '\n' to move to a new line.

clear all;
close all;
clc;

filepath_metadata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\Archive\metadata.xlsx';

% Define case to simulate
case2sim = '180817_1';

% Pressure loss array
dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04];

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata);
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
cases = readtable(filepath_metadata,opts);

% Select 
case2sim = cases(strcmp(cases.ID,case2sim),:);

% Estimate mdot based on dpdx (Blasius)
vel = (dpdx.^4.*case2sim.d_h.^5./(1./16.*0.316.^4.*case2sim.K_PL_scaled.*case2sim.rho_scaled.^3)).^(1./7);
A = pi./4.*(case2sim.d_o.^2-case2sim.d_i.^2);
Q = vel.*A;
mdot = Q.*case2sim.rho_scaled;


% Parameters
if strcmp(case2sim.Geometry,'Pipe')
    geo = 'p'; 
else
    geo = 'a';
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

% Path and filename of journalfile
filepath_journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\remote';
fileID = fopen([filepath_journal '\RUN_PressureGradientSweep.jou'],'w');


%% Write journal file
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '; Pressure gradient sweep of turbulent pipe/annular flow\r\n');
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '; 14.06.2018, alexander.busch@alumni.ntnu.no\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '; Reads generic case and sweeps pressure gradient\r\n');
fprintf(fileID, '; Attention: Change absolute path in rename command to suit working directory \r\n');
fprintf(fileID, '; Paste "file/read-journal RUN_PressureGradientSweep.jou" to run simulations for defined pressure gradients\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '; Read case and change color scheme\r\n');
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
if strcmp(geo,'p')
    fprintf(fileID, '/file/read-case-data pipe.cas OK\r\n');
else
    fprintf(fileID, '/file/read-case-data annulus.cas OK\r\n');    
end
fprintf(fileID, '/display/set/colors/color-scheme classic\r\n');
fprintf(fileID, '/display/set/colors/background "white"\r\n');
fprintf(fileID, '/display/set/colors/foreground "black"\r\n');
fprintf(fileID, '/display/re-render\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '; Purge working directory\r\n');
fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '!if exist reports rmdir /s /q reports\r\n');
fprintf(fileID, '!mkdir reports\r\n');
fprintf(fileID, '!if exist exports rmdir /s /q exports\r\n');
fprintf(fileID, '!mkdir exports\r\n');
fprintf(fileID, '!if exist monitors rmdir /s /q monitors\r\n');
fprintf(fileID, '!mkdir monitors\r\n');
fprintf(fileID, '!if exist solved rmdir /s /q solved\r\n');
fprintf(fileID, '!mkdir solved\r\n');
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

% Loop all dpdx
for ii = 1:length(dpdx)
    if strcmp(turbmod,'epsilon')
        fprintf(fileID, ['/define/periodic-conditions/massflow-rate-specification ' num2str(mdot(ii)) ' -' num2str(dpdx(ii)) ' 0.5 5 1 0 0\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 3 y y "monitors/' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '_mass-flow-rate.txt" 1\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 3 y y "monitors/' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '_bulk-velocity.txt" 1\r\n']);
        fprintf(fileID, '/solve/monitors/statistic/monitors y\r\n'); 
        fprintf(fileID, '/solve/monitors/statistic/print y\r\n'); 
        fprintf(fileID, '/solve/monitors/statistic/plot y\r\n'); 
        fprintf(fileID, '/solve/monitors/statistic/write y\r\n'); 
        fprintf(fileID, '/solve/monitors/statistic/window 4\r\n');
    else
        fprintf(fileID, ['/define/periodic-conditions/pressure-gradient-specification -' num2str(dpdx(ii)) ' 1 0 0\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 3 y y "monitors/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '_mass-flow-rate.txt" 1\r\n']);
        fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 3 y y "monitors/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '_bulk-velocity.txt" 1\r\n']);
    end
    fprintf(fileID, ['/solve/initialize/set-defaults/x-velocity ' num2str(vel(ii)) '\r\n']);
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
    fprintf(fileID, '/solve/iterate 20000\r\n');
%     fprintf(fileID, ['/file/export/ascii "exports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate specific-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
%     fprintf(fileID, ['/report/summary y "reports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt"\r\n']);
%     fprintf(fileID, ['/file/write-case-data "solved/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.cas"']);
    if strcmp(turbmod,'epsilon')
        fprintf(fileID, ['/file/export/ascii "exports/' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
        fprintf(fileID, ['/report/summary y "reports/' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '.txt"\r\n']);
        fprintf(fileID, ['/file/write-case-data "solved/' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '.cas"\r\n']);
        fprintf(fileID, '!move /y statistics-per-pr-grad.out monitors\\statistics-per-pr-grad.out\r\n');
        fprintf(fileID, ['!rename ' strrep(filepath_journal,'\','\\') '\\monitors\\statistics-per-pr-grad.out "' geo '_' fluid '-scaled_mdot=' num2str(mdot(ii)) '_dpdx.txt"\r\n']);
    else
        fprintf(fileID, ['/file/export/ascii "exports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt" line_xy () , velocity-magnitude x-velocity y-velocity viscosity-turb turb-kinetic-energy turb-diss-rate specific-diss-rate turb-reynolds-number-rey production-of-k y-plus () ,\r\n']);
        fprintf(fileID, ['/report/summary y "reports/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.txt"\r\n']);
        fprintf(fileID, ['/file/write-case-data "solved/' geo '_' fluid '-scaled_dpdx=' num2str(dpdx(ii)) '.cas"\r\n']);
    end
    fprintf(fileID, '\r\n');
end

fclose(fileID);