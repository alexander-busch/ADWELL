% Orders of magnitude estimate of non-dimensional numbers
% characterizing cuttings transport

Initialize;
Parameters;

%% Rheological properties/flowcurve at surface

% Compute Fann dial readings
DR = (YP+PV/(1022-511)*SR_Fann); % [°]

% Pipe shear rate range PL coefficients (API RP 13D)
% --> PL fluid behavior for both pipe and annular flow
n_P = 3.32*log10((DR(:,1))./DR(:,2)); n_P(1:4) = 1;
K_P = DR(:,2)./(511.^n_P);

% Annulus shear rate range PL coefficients (API RP 13D),
% evaluated based on YP/PV fit with above computed Fann DR for low SR range
% ---> Bingham fluid behavior for annular flow
n_A = 0.657*log10(DR(:,4)./DR(:,6)); n_A(1:4) = 1;
K_A = DR(:,4)./(170.3.^n_A);


%% Rheological properties/flowcurve downhole



%% Bulk velocities

U = [Q(:,1)./A(:,1)...  % Min. ax. pipe velocity
    Q(:,2)./A(:,1)...   % Max. ax. pipe velocity
    Q(:,1)./A(:,2)...   % Min. ax. annular velocity, 0 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,1)/60).^2).^0.5...   % Min. ax. annular velocity, 1 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,2)/60).^2).^0.5...	% Min. ax. annular velocity, 2 rpm
    Q(:,2)./A(:,2)...   % Max. ax. annular velocity, 0 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,1)/60).^2).^0.5...   % Max. ax. annular velocity, 1 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,2)/60).^2).^0.5];	% M4ax. ax. annular velocity, 2 rpm


%% 

% Newtonian wall shear rate 
SR_N = [8*U(:,1:2)./d_h(:,1) 12*U(:,3:8)./d_h(:,2)];

% Pre-calc
Re_nom = rho_f*[U(:,1:2).*d_h(:,1) U(:,3:8).*d_h(:,2)];
Sh_denom = ((rho_s-rho_f)*9.81*d_p(2));
   
% Initialize variables
SR_PL = zeros(length(d_h), length(d_h), length(PV));
SR_YPPV = zeros(length(d_h), length(d_h), length(PV));
eta_f_PL = zeros(length(d_h), length(d_h), length(PV));
eta_f_YPPV = zeros(length(d_h), length(d_h), length(PV));
Re_PL = zeros(length(d_h), length(d_h), length(PV));
Re_YPPV = zeros(length(d_h), length(d_h), length(PV));
Sh_PL = zeros(length(d_h), length(d_h), length(PV));
Sh_YPPV = zeros(length(d_h), length(d_h), length(PV));
% Loop all fluids
for fluid=1:length(PV)
    % True wall shear rate based on local n' and K' of Metzner & Reed as
    % given in API RP 13D for pipe (higher shear rate range --> PL = YPPV)
    % and annulus (lower shear rate range --> PL not equal YPPV).
    SR_PL(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid))) * SR_N(:,1:2),...
        ((1+2*n_P(fluid))/(3*n_P(fluid))) * SR_N(:,3:8)];
    SR_YPPV(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid))) * SR_N(:,1:2),...
        ((1+2*n_A(fluid))/(3*n_A(fluid))) * SR_N(:,3:8)];

    % PL model evaluated with true shear rate and high (pipe) shear
    % rate range coefficients (--> Ostwald fluid)
    eta_f_PL(:,:,fluid) = lbf100f22Pa.*FannDEG2lbf100f2.*...
        (K_P(fluid).*SR_PL(:,:,fluid).^n_P(fluid))./SR_PL(:,:,fluid);
    % PL model evaluated with true shear rate and high and low (pipe and
    % annular) shear rate range (--> Bingham fluid)
    eta_f_YPPV(:,:,fluid) = lbf100f22Pa.*FannDEG2lbf100f2.*...
        (K_A(fluid).*SR_YPPV(:,:,fluid).^n_A(fluid))./SR_YPPV(:,:,fluid);
end

% Re number
Re_PL = Re_nom./eta_f_PL;
Re_YPPV = Re_nom./eta_f_YPPV;

% Sh number
Sh_PL = eta_f_PL.*SR_PL./Sh_denom;
Sh_YPPV = eta_f_YPPV.*SR_YPPV./Sh_denom;


%% Old computation with different shear rate estimates

% Select fluid 
fluid = 8;
% Newtonian wall shear rate 
SR_N = [8*U(:,1:2)./d_h(:,1)...
    12*U(:,3:8)./d_h(:,2)];
% True wall shear rate based on local n' and K' of Metzner & Reed as given
% in API RP 13D for pipe (higher shear rate range --> PL = YPPV) and
% annulus (lower shear rate range --> PL not equal YPPV).
SR_P =      ((1+3*n_P(fluid))/(4*n_P(fluid))) * SR_N(:,1:2);
SR_A_PL =   ((1+2*n_P(fluid))/(3*n_P(fluid))) * SR_N(:,3:8);
SR_A_YPPV = ((1+2*n_A(fluid))/(3*n_A(fluid))) * SR_N(:,3:8);

% eta_f is a three-dimensional matrix, where the
% - first dimension are the eight hole sections (columns)
% - second dimension are the min and max values for pipe, annulus, rotation
% - third dimension corresponds to the following list.
eta_f = zeros(8,8,7);

% PL model evaluated with Newtonian shear rate and high (pipe) shear
% rate range coefficients
eta_f(:,:,1) = lbf100f22Pa.*FannDEG2lbf100f2.* (K_P(fluid).*SR_N.^n_P(fluid))./SR_N;
% PL model evaluated with Newtonian shear rate and high and low (pipe and
% annular) shear rate range coefficients
eta_f(:,:,2) = lbf100f22Pa.*FannDEG2lbf100f2.* (K_A(fluid).*SR_N.^n_A(fluid))./SR_N;
% PL model evaluated with true shear rate and high (pipe) shear
% rate range coefficients (--> Ostwald fluid)
eta_f(:,:,3) = lbf100f22Pa.*FannDEG2lbf100f2.* (K_P(fluid).*[SR_P SR_A_PL].^n_P(fluid))./[SR_P SR_A_PL];
% PL model evaluated with true shear rate and high and low (pipe and
% annular) shear rate range (--> Bingham fluid)
eta_f(:,:,4) = lbf100f22Pa.*FannDEG2lbf100f2.* (K_A(fluid).*[SR_P SR_A_YPPV].^n_A(fluid))./[SR_P SR_A_YPPV];

% YPPV model evaluated with Newtonian shear rate and high (pipe) shear
% rate range coefficients
eta_f(:,:,5) = lbf100f22Pa.*FannDEG2lbf100f2.*(YP(fluid).*ones(size(SR_N))+PV(fluid)./(1022-511).*SR_N)./SR_N;
% YPPV model evaluated with true shear rate and high (pipe) shear
% rate range coefficients (--> Ostwald fluid)
eta_f(:,:,6) = lbf100f22Pa.*FannDEG2lbf100f2.*(YP(fluid).*ones(size([SR_P SR_A_PL]))+PV(fluid)./(1022-511).*[SR_P SR_A_PL])./[SR_P SR_A_PL];
% YPPV model evaluated with true shear rate and high and low (pipe and
% annular) shear rate range (--> Bingham fluid)
eta_f(:,:,7) = lbf100f22Pa.*FannDEG2lbf100f2.*(YP(fluid).*ones(size([SR_P SR_A_YPPV]))+PV(fluid)./(1022-511).*[SR_P SR_A_YPPV])./[SR_P SR_A_YPPV];

% Ostwald vs. Bingham behavior based on PL model
% OvsB_PL = eta_f(:,:,3)./eta_f(:,:,4);

% Ostwald vs. Bingham behavior based on YPPV model
% OvsB_YPPV = eta_f(:,:,6)./eta_f(:,:,7);

% PL vs. YPPV for Ostwald behavior
% PLvsYPPV_O = eta_f(:,:,3)./eta_f(:,:,6);

% PL vs. YPPV for Bingham behavior
% PLvsYPPV_B = eta_f(:,:,4)./eta_f(:,:,7);

% Re is a three-dimensional matrix where the third dimension corresponds to
% the list of viscosities as given above.
Re = rho_f*[U(:,1:2).*d_h(:,1) U(:,3:8).*d_h(:,2)]./eta_f;

% Ostwald vs. Bingham behavior based on PL model
% OvsB_PL = Re(:,:,2)./Re(:,:,3);

% Ostwald vs. Bingham behavior based on YPPV model
% OvsB_YPPV = Re(:,:,5)./Re(:,:,6);

% PL vs. YPPV for Ostwald behavior
% PLvsYPPV_O = Re(:,:,2)./Re(:,:,5);

% PL vs. YPPV for Bingham behavior
% PLvsYPPV_B = Re(:,:,3)./Re(:,:,6);

% PL model evaluated with true shear rate and high (pipe) shear
% rate range coefficients (--> Ostwald fluid)
Sh = eta_f(:,:,3).*[SR_P SR_A_PL]./((rho_s-rho_f)*9.81*d_p(2));
Sh = eta_f(:,:,4).*[SR_P SR_A_YPPV]./((rho_s-rho_f)*9.81*d_p(2));


%% Particle settling

% Iterate settling velocity and apparent viscosity
Re_p = rho_f.*U.*d_h./eta_f;



%% Plot

PlotParameters;
PlotFlowCurves;
PlotBulkVelocities;
PlotShearRates;
PlotApparentViscosities;
PlotReynoldsNumbers;
PlotShieldsNumbers;