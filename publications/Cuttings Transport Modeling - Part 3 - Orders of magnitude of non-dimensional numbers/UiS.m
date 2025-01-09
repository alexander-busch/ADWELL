%% Geometry [in]
% First row is the inner drill pipe diameter,
% second row is the bit/hole diameter and the outer drill pipe diameter
d_o = [0.040; 0.040]';
d_i = [0.000; 0.025]';
% Hydraulic diameter
d_h = d_o-d_i;
% Cross-sectional area
A = pi/4*(d_o.^2-d_i.^2);


%% Flow rates [gpm]
U_sl = [0.1 1.5];
Q = U_sl*A(1);


%% Drill string rotation rates [rpm]
rpm = [70 160];



%% Material properties

% Particle diameter
d_p = [0.3 1.2]/1000;

% Densities
rho_f = 1000;
rho_s = 2500;

% Drilling fluid rheological properties
mu_inf = 0.00102; % [Pa.s] H2O
mu_0 = 0.026; % [Pa.s] H2O
lambda_Cr = 0.008;
n_Cr = 0.37;


%% No cuttings bed present

U = [Q(:,1)./A(:,1)...  % Min. ax. pipe velocity
    Q(:,2)./A(:,1)...   % Max. ax. pipe velocity
    Q(:,1)./A(:,2)...   % Min. ax. annular velocity, 0 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,1)/60).^2).^0.5...   % Min. ax. annular velocity, 1 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,2)/60).^2).^0.5...	% Min. ax. annular velocity, 2 rpm
    Q(:,2)./A(:,2)...   % Max. ax. annular velocity, 0 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,1)/60).^2).^0.5...   % Max. ax. annular velocity, 1 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_h(:,2)/2.*rpm(:,2)/60).^2).^0.5];	% M4ax. ax. annular velocity, 2 rpm

% Newtonian wall shear rate 
SR_N = [8*U(:,1:2)./d_h(:,1) 12*U(:,3:8)./d_h(:,2)];

% Pre-calc
Re_nom = rho_f*[U(:,1:2).*d_h(:,1) U(:,3:8).*d_h(:,2)];
Sh_denom = ((rho_s-rho_f)*9.81*d_p);
Ar_nom = rho_f*(rho_s-rho_f)*9.81*d_p.^3;

for fluid=1:2
    
    if fluid == 1
        n_PL = ones(1,length(SR_N));
        K_PL = mu_inf*ones(1,length(SR_N));
    else
        addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';
        [ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR_N );
    end
    
        
    
    
    
    
    
    % True wall shear rate based on local n' and K' of Metzner & Reed as
    % given in API RP 13D for pipe (higher shear rate range --> PL = YPPV)
    % and annulus (lower shear rate range --> PL not equal YPPV).
    SR_PL(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid))) * SR_N(:,1:2),...
        ((1+2*n_P(fluid))/(3*n_P(fluid))) * SR_N(:,3:8)];
    SR_YPPV(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid))) * SR_N(:,1:2),...
        ((1+2*n_A(fluid))/(3*n_A(fluid))) * SR_N(:,3:8)];

    % PL model evaluated with true shear rate and high (pipe) shear
    % rate range coefficients (--> Ostwald fluid)
    eta_f_PL(:,:,fluid) = (K_P(fluid).*SR_PL(:,:,fluid).^n_P(fluid))./SR_PL(:,:,fluid);
    % PL model evaluated with true shear rate and high and low (pipe and
    % annular) shear rate range (--> Bingham fluid)
    eta_f_YPPV(:,:,fluid) = (K_A(fluid).*SR_YPPV(:,:,fluid).^n_A(fluid))./SR_YPPV(:,:,fluid);
    % Bingham viscosity
    my(:,:,fluid) = PV(fluid);
    
    Yi(:,:,fluid) = YP(fluid).*[d_h(:,1).*ones(length(d_h),2) d_h(:,2).*ones(length(d_h),6)]./(PV(fluid).*U);
    
end