% Read metadata

clear all;
close all;
clc;

filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive';
filepath_metadata = [filepath_simdata  '\' 'metadata.xlsx'];

% Import table with units and description fields
opts = detectImportOptions(filepath_metadata);
opts.VariableUnitsRange='A2';
opts.VariableDescriptionsRange='A3';
opts = setvartype(opts,{'dpdx', 'mdot_f', 'Q_f', 'U_f_x', 'rpm'},'char');
% disp([opts.VariableNames' opts.VariableTypes']);
cases = readtable(filepath_metadata,opts);

% Conversion from character strings to numerical arrays
for aa=1:height(cases)
    cases.dpdx{aa} = str2num(cases.dpdx{aa});
    cases.mdot_f{aa} = str2num(cases.mdot_f{aa});
    cases.Q_f{aa} = str2num(cases.Q_f{aa});
    cases.U_f_x{aa} = str2num(cases.U_f_x{aa});
    cases.rpm{aa} = str2num(cases.rpm{aa});
end
    dpdx = cases.dpdx{aa};


% Loop all table rows
for aa=1:height(cases)
    % aa = 7;
 

    % Geometry coefficient beta
    if strcmp(cases.Geometry(aa),'Pipe')
        cases.beta(aa) = 3;
    else
        cases.beta(aa) = 2;
    end   

    % Newtonian shear rate
    cases.SR_newt{aa,1} = 24./cases.beta(aa).*cases.U_f_x{aa,1}./cases.d_h(aa);
    
    % PL shear rate, handle n = 1 seperately due to formulation of PL shear
    % rate correction where for n_PL = 1 the denominator in the exponent
    % becomes zero
    if (cases.n_PL_f(aa)==1)
        cases.SR_PL{aa,1} = cases.SR_newt{aa,1};
    else
        cases.SR_PL{aa,1} = cases.SR_newt{aa,1} .* ((3.*cases.n_PL_f(aa)+1)./(4.*cases.n_PL_f(aa))) .^ (cases.n_PL_f(aa)./(cases.n_PL_f(aa)-1));
    end
    
    % PL viscosity
    cases.eta_PL{aa,1} = cases.K_PL_f(aa).*cases.SR_PL{aa,1}.^(cases.n_PL_f(aa)-1);
    
    % Reynolds number based on apparent viscosity
    % cases.Re_G{aa,1} = cases.rho_f(aa).*cases.U_f_x{aa,1}.*cases.d_h(aa)./cases.eta_PL{aa,1};
    
    % Generalized Reynolds number (Metzner-Reed 1955)
    % cases.Re_MR{aa,1} = cases.rho_f(aa).*cases.U_f_x{aa,1}.^(2-cases.n_PL_f(aa)).*cases.d_h(aa).^cases.n_PL_f(aa)./(cases.K_PL_f(aa).* ((3.*cases.n_PL_f(aa)+1)./(4.*cases.n_PL_f(aa))).^cases.n_PL_f(aa) .*8.^(cases.n_PL_f(aa)-1));
    
    % Generalized Reynolds number (Metzner-Reed 1955) with annular geometry correction of Kozicki et al. (1966) and Delplace and Leuliet (1995) 
    cases.Re_G{aa,1} = cases.rho_f(aa).*cases.U_f_x{aa,1}.^(2-cases.n_PL_f(aa)).*cases.d_h(aa).^cases.n_PL_f(aa)./((24./cases.beta(aa)).^(cases.n_PL_f(aa)-1).*cases.K_PL_f(aa).*((1+cases.n_PL_f(aa).*cases.beta(aa))./(cases.n_PL_f(aa)+cases.n_PL_f(aa).*cases.beta(aa))).^cases.n_PL_f(aa));
    
     % Wall shear stress
    cases.tauw{aa,1} = cases.dpdx{aa,1}.*cases.d_h(aa)./4;

    % Dynamic pressure
    cases.p_dyn{aa,1} = cases.rho_f(aa).*cases.U_f_x{aa,1}.^2./2;
    
    % Friction factor
    cases.f{aa,1} = cases.tauw{aa,1}./cases.p_dyn{aa,1};

    % Friction velocity
    cases.ufricf{aa,1} = sqrt(cases.tauw{aa,1}./cases.rho_f(aa));

    % Wall shear rate
    cases.SR_wall{aa,1} = (cases.tauw{aa,1}./cases.K_PL_f(aa)).^(1./cases.n_PL_f(aa));

    % Wall viscosity
    cases.eta_wall{aa,1} = cases.tauw{aa,1}./cases.SR_wall{aa,1};


    for bb=1:length(cases.ufricf{aa,1})

        % u+
       %cases.u_plus{aa,1}(bb,:)= cases.u_x{aa,1}(bb,:)./cases.ufricf{aa,1}(bb);

        % y+
        %cases.y_plus{aa,1}(bb,:) = (cases.d_o(aa)/2 - cases.y{aa,1}(bb,:)) .* cases.rho_f(aa).*cases.ufricf{aa,1}(bb)./cases.eta_wall{aa,1}(bb);

    end
    
    % Courant number 
    % cases.dt

end

save('wellboredata.mat', 'cases');



