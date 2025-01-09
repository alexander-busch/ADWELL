% Create table of experimental Flow Loop data and Pipe Viscometry variables

%% Create Pipe Viscometry  table with subset (Date, Fluid, Usl_p, dpdl_p) of Exp
PipeVis = Exp(:,[1 2 10 8]);
 

%% Compute Newtonian Shear Rate and Shear stress at the Wall
PipeVis.NewtSR = NewtShearRate(do, PipeVis.Usl_p);
PipeVis.WallSS = WallShearStress(do, PipeVis.dpdl_p);


%% Initialize n & K
PipeVis.KPrime = zeros(height(PipeVis),1);
PipeVis.nPrime = zeros(height(PipeVis),1);

PipeVis.K = zeros(height(PipeVis),1);
PipeVis.n = zeros(height(PipeVis),1);

%% Compute n' & K'  for each fluid
for i=1:length(fluidlist) % Loop through fluid cases
    fluidname = fluidlist{i}; % Define casename from caselist (Cell element)
    filter = strcmp(PipeVis.Fluid,fluidname);
    nPrime=zeros(nnz(filter),1);
    KPrime=zeros(nnz(filter),1);
        
    if strcmp(fluidname, 'H2O') % Assign H2O standard values

        
        nPrime(:,1) = 1;
        KPrime(:,1) = 0.00102;

        PipeVis.nPrime(filter) = nPrime;
        PipeVis.KPrime(filter) = KPrime;

    else

        if strcmp(fluidname, 'PAC2') % Redefine filter to account for lates PAC2 data
            filter = strcmp(PipeVis.Fluid,fluidname) & datetime(PipeVis.Date)=='29-Jul-2016'; % strcmp(PipeVis.Date,'2016-07-29');
            nPrime=zeros(nnz(filter),1);
            KPrime=zeros(nnz(filter),1);
        end

        % Create vectors used for curve fitting
        NeSR = PipeVis.NewtSR(filter);
        TauW = PipeVis.WallSS(filter);

        if strcmp(fluidname, 'PAC1') % Create laminar dataset for NeSR and TauW

            LasDatPoi = length(NeSR) - 5; % Create lam. dataset by substracting last data points which are turbulent acc. to f-ReG-plot
            NeSR = NeSR(1:LasDatPoi);
            TauW = TauW(1:LasDatPoi);

        end 

        fittedmodel = CreatePLFit(NeSR, TauW);
        coefficients = coeffvalues(fittedmodel);

        nPrime(:,1) = coefficients(2);
        KPrime(:,1) = coefficients(1);

        PipeVis.nPrime(filter) = nPrime;
        PipeVis.KPrime(filter) = KPrime;
    
    end
    
end
    
%% Compute n' & K'
PipeVis.n = PipeVis.nPrime;
PipeVis.K = PipeVis.KPrime./(((3.*PipeVis.n+1)./(4.*PipeVis.n)).^PipeVis.n);


%% Clean up
clear nPrime KPrime coefficients fittedmodel LasDatPoi NeSR TauW