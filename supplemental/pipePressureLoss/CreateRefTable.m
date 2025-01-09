% Reference values - Analytical and empirical solutions


%% Number of elements for each fluid case
elements = 51;


%% Initialize vectors
Fluid = cell(length(fluidlist)*elements,1);
MassFlowRate = zeros(length(fluidlist)*elements,1);
RPM = zeros(length(fluidlist)*elements,1);
Density = zeros(length(fluidlist)*elements,1);
K = zeros(length(fluidlist)*elements,1);
n = zeros(length(fluidlist)*elements,1);


%% Fill up vectors with respective values
for i = 1:length(fluidlist)
    fluidname = fluidlist(i); % Get current fluid name
    filter = strcmp(Exp.Fluid,fluidname); % Define filter for Exp table
    
    % Fill up previously defined vectors with respective values
    
    Fluid(i*elements-elements+1:i*elements) = fluidname;
    
    MassFlowRate(i*elements-elements+1:i*elements) = linspace(0,2,elements);
  
    RPM;
    
    Density(i*elements-elements+1:i*elements) = mean(Exp.Density(filter));

    K(i*elements-elements+1:i*elements) = mean(Exp.K(filter));
    n(i*elements-elements+1:i*elements) = mean(Exp.n(filter));
  
end


%% Create table & cleanup
Ref = table(Fluid, MassFlowRate, RPM, Density, K, n); % Create ref table with fluid and massflowrate
clear fluidname Fluid_i Fluid filter MassFlowRate Density K n;


%% Compute derived variables
Ref.Usl_p = SuperficialVelocity (Ref.MassFlowRate, Ref.Density, do, di, 3);
Ref.Usl_a = SuperficialVelocity (Ref.MassFlowRate, Ref.Density, do, di, 2);
Ref.ReG_p = GeneralizedReynoldsNumber ( Ref.Density, Ref.Usl_p, do, di, 3, Ref.n, Ref.K );
Ref.ReG_a = GeneralizedReynoldsNumber ( Ref.Density, Ref.Usl_p, do, di, 2, Ref.n, Ref.K );
Ref.f_p = FanningFrictionFactor_Ref ( Ref.ReG_p, Ref.Fluid, 3 );
Ref.f_a = FanningFrictionFactor_Ref ( Ref.ReG_a, Ref.Fluid, 2 );
% Ref.f_ae = FanningFrictionFactor_Ref  ( Exp.ReG_a, Ref.Fluid, 2 );
Ref.dpdl_p = PressureGradient ( Ref.f_p, Ref.Density, Ref.Usl_p, do, di, 3);
Ref.dpdl_a = PressureGradient ( Ref.f_a, Ref.Density, Ref.Usl_a, do, di, 2);
% Ref.dpdl_ae = PressureGradient ( Ref.f_ae, Exp.Density, Exp.Usl_ae, do, di, 2);