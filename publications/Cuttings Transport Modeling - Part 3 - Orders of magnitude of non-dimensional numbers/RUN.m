% Orders of magnitude estimate of non-dimensional numbers
% characterizing cuttings transport

Init;
Parameters;

fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Turbulence level estimate\figures\';


%% Rheological properties/flowcurve at surface based on API RP 13D

% Shear stress [Pa] as function of shear rate [1/s]
tau_Fann = YP+PV.*SR_Fann; %
% plot(SR_Fann,tau)

% Pipe shear rate range PL coefficients (API RP 13D)
% --> PL fluid behavior for both pipe and annular flow
n_P = 3.32*log10((tau_Fann(:,1))./tau_Fann(:,2)); n_P(1:4) = 1; % [-]
K_P = tau_Fann(:,2)./(511.^n_P); % [Pa.s^n]

% Annulus shear rate range PL coefficients (API RP 13D),
% evaluated based on YP/PV fit with above computed Fann DR for low SR range
% ---> Bingham fluid behavior for annular flow
n_A = 0.657*log10(tau_Fann(:,4)./tau_Fann(:,6)); n_A(1:4) = 1; % [-]
K_A = tau_Fann(:,4)./(170.3.^n_A); % [Pa.s^n]



%% Rheological properties/flowcurve downhole



%% No cuttings bed present

U = [Q(:,1)./A(:,1)...  % Min. ax. pipe velocity
    Q(:,2)./A(:,1)...   % Max. ax. pipe velocity
    Q(:,1)./A(:,2)...   % Min. ax. annular velocity, 0 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_i(:,2).*rpm(:,1)/60).^2).^0.5...   % Min. ax. annular velocity, 1 rpm
    ((Q(:,1)./A(:,2)).^2+(0.5*2*pi*d_i(:,2).*rpm(:,2)/60).^2).^0.5...	% Min. ax. annular velocity, 2 rpm
    Q(:,2)./A(:,2)...   % Max. ax. annular velocity, 0 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_i(:,2).*rpm(:,1)/60).^2).^0.5...   % Max. ax. annular velocity, 1 rpm
    ((Q(:,2)./A(:,2)).^2+(0.5*2*pi*d_i(:,2).*rpm(:,2)/60).^2).^0.5];	% M4ax. ax. annular velocity, 2 rpm

% Newtonian wall shear rate 
SR_N = [8*U(:,1:2)./d_h(:,1) 12*U(:,3:8)./d_h(:,2)];

% Pre-calc
Re_nom = rho_f*[U(:,1:2).*d_h(:,1) U(:,3:8).*d_h(:,2)];
Sh_denom = ((rho_s-rho_f)*9.81*d_p);
Ar_nom = rho_f*(rho_s-rho_f)*9.81*d_p.^3;
   
% Initialize variables
columns=size(U);
SR_PL = zeros(length(d_h) , columns(2), length(PV));
SR_YPPV = zeros(length(d_h) , columns(2), length(PV));
eta_f_PL = zeros(length(d_h) , columns(2), length(PV));
eta_f_YPPV = zeros(length(d_h) , columns(2), length(PV));
my = zeros(length(d_h) , columns(2), length(PV));
Re_PL = zeros(length(d_h) , columns(2), length(PV));
Re_YPPV = zeros(length(d_h) , columns(2), length(PV));
Re_Bi = zeros(length(d_h) , columns(2), length(PV));
Yi = zeros(length(d_h) , columns(2), length(PV));
Sh_PL = zeros(length(d_h) , columns(2), length(PV));
Sh_YPPV = zeros(length(d_h) , columns(2), length(PV));
% Loop all fluids
for fluid=1:length(PV)
    % True wall shear rate based on local n' and K' of Metzner & Reed as
    % given in API RP 13D for pipe (higher shear rate range --> PL = YPPV)
    % and annulus (lower shear rate range --> PL not equal YPPV).
    
    if n_P(fluid)==1
        SR_PL(:,:,fluid) = [1 * SR_N(:,1:2),...
            1 * SR_N(:,3:8)];
        SR_YPPV(:,:,fluid) = [1 * SR_N(:,1:2),...
            1 * SR_N(:,3:8)];
    else
        
        SR_PL(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,1:2),...
            ((1+2*n_P(fluid))/(3*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,3:8)];
        SR_YPPV(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,1:2),...
            ((1+2*n_A(fluid))/(3*n_A(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,3:8)];
    end
    
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

% Re number
Re_PL = Re_nom./eta_f_PL; Re_YPPV = Re_nom./eta_f_YPPV;
Re_cr_PL = [3250-1150*n_P 3250-1150*n_A];

Re_Bi = Re_nom./my;
He = Re_Bi.*Yi;

% Compute critical Reynolds number element-wise
Re_cr_Bi = zeros(size(He));
for fluid = 1:length(PV)
    for jj=1:min(size(He))
        for kk=1:min(size(He))
            % Alternative: Hanks (1963) as given in Darby and Melson (1981) or
            % Chhabra and Richardson (2008)
            if He(jj,kk,fluid)==0
                Re_cr_Bi(jj,kk,fluid)=2100;
            elseif He(jj,kk,fluid)<1700
                Re_cr_Bi(jj,kk,fluid)=2100./(1+8.3e-8*log10(He(jj,kk,fluid)).^13);
            elseif He(jj,kk,fluid)<5e4
                Re_cr_Bi(jj,kk,fluid)=80.*He(jj,kk,fluid).^0.4;
            else
                Re_cr_Bi(jj,kk,fluid)=25.*He(jj,kk,fluid).^0.5;
            end
        end
    end
end


% Select fluid and particle example
fluid = 8;
particle = 2;


%Particle Re number
Re_p = rho_f*U.*d_p(particle)./eta_f_PL; Re_p(:,1:2,:)=0;

% Critical Shields number (Soulsby 1997)
d_star = (9.81.*(rho_s./rho_f-1)./(eta_f_PL./rho_f).^2).^(1/3).*d_p(particle); d_star(:,1:2,:)=0;
Theta_cr = 0.24./d_star+0.055.*(1-exp(-0.02.*d_star)); Theta_cr(:,1:2,:)=0;

% Sh number
Sh_PL = eta_f_PL.*SR_PL./Sh_denom(particle); Sh_PL(:,1:2,:)=0;
Sh_YPPV = eta_f_YPPV.*SR_YPPV./Sh_denom(particle); Sh_YPPV(:,1:2,:)=0;

% Normalize Sh number w/ critical Shields number
Sh_PL = Sh_PL./Theta_cr;
Sh_YPPV = Sh_YPPV./Theta_cr;

% Plot

PlotParameters;
PlotFlowCurves;
PlotBulkVelocities;
PlotShearRates;
PlotApparentViscosities;
PlotReynoldsNumbers;
PlotShieldsNumbers;





% Remove everything larger than 26" for non-vertical analysis as recomended
% by TAG
d_i(1:3,:) = [];
d_o(1:3,:) = [];
d_h(1:3,:) = [];
A(1:3,:) = [];
Q(1:3,:) = [];
rpm(1:3,:) = [];
%% Cuttings bed present

fig_Re_Sh = CreateFigure( '', '', '', {'log','log'}, 'DINA4' );
axis off;
fig_Re_Sh_sub = tight_subplot(3,3,[.1 .06],[.06 .03],[.06 .02]);

fig_Ouriemi = CreateFigure( '', '', '', {'log','log'}, 'DINA4' );
axis off;
fig_Ouriemi_sub = tight_subplot(3,3,[.1 .06],[.06 .03],[.06 .02]);



counter=1;
% Bed height normalized with annular outer radius
% h_b = [0.00 0.33 0.66];
h_b = [0.25 0.50 0.75];





% Loop all bed heights
% ii = 2;
for ii = 1:length(h_b)

	% Compute effective variables for annulus and pipe for current bed
    % height
    if h_b(ii)>0
        [ w_a, w_p, h_f, A_a_f, A_p_f, d_h_a_eff, d_h_p_eff ] = CompEffFlowVar( d_i(:,2), d_o(:,2), d_h(:,2), h_b(ii)*d_o(:,2)/2 );
        
        d_h = [d_h(:,1) d_h_a_eff];
        A =  [A(:,1) A_a_f];
    else

        h_f = 0.001*ones(length(d_h(:,2)),1);

    end 
    
    
    % Bulk velocities
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

    % Initialize variables
    SR_PL = zeros(length(d_h) , columns(2), length(PV));
    SR_YPPV = zeros(length(d_h) , columns(2), length(PV));
    eta_f_PL = zeros(length(d_h) , columns(2), length(PV));
    eta_f_YPPV = zeros(length(d_h) , columns(2), length(PV));
    my = zeros(length(d_h) , columns(2), length(PV));
    Re_PL = zeros(length(d_h) , columns(2), length(PV));
    Re_YPPV = zeros(length(d_h) , columns(2), length(PV));
    Re_Bi = zeros(length(d_h) , columns(2), length(PV));
    Yi = zeros(length(d_h) , columns(2), length(PV));
    Sh_PL = zeros(length(d_h) , columns(2), length(PV));
    Sh_YPPV = zeros(length(d_h) , columns(2), length(PV));
    % Loop all fluids
    for fluid=1:length(PV)
        % True wall shear rate based on local n' and K' of Metzner & Reed as
        % given in API RP 13D for pipe (higher shear rate range --> PL = YPPV)
        % and annulus (lower shear rate range --> PL not equal YPPV).
        
        if n_P(fluid)==1
            SR_PL(:,:,fluid) = [1 * SR_N(:,1:2),...
                1 * SR_N(:,3:8)];
            SR_YPPV(:,:,fluid) = [1 * SR_N(:,1:2),...
                1 * SR_N(:,3:8)];
        else

            SR_PL(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,1:2),...
                ((1+2*n_P(fluid))/(3*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,3:8)];
            SR_YPPV(:,:,fluid) = [((1+3*n_P(fluid))/(4*n_P(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,1:2),...
                ((1+2*n_A(fluid))/(3*n_A(fluid)))^(n_P(fluid)/(n_P(fluid)-1)) * SR_N(:,3:8)];
        end
            
        % PL model evaluated with true shear rate and high (pipe) shear
        % rate range coefficients (--> Ostwald fluid)
        eta_f_PL(:,:,fluid) = (K_P(fluid).*SR_PL(:,:,fluid).^n_P(fluid))./SR_PL(:,:,fluid);
        % PL model evaluated with true shear rate and high and low (pipe and
        % annular) shear rate range (--> Bingham fluid)
        eta_f_YPPV(:,:,fluid) = (K_A(fluid).*SR_YPPV(:,:,fluid).^n_A(fluid))./SR_YPPV(:,:,fluid);
        % Bingham viscosity
        my(:,:,fluid) = PV(fluid);
        % Yield number
        Yi(:,:,fluid) = YP(fluid).*[d_h(:,1).*ones(length(d_h),2) d_h(:,2).*ones(length(d_h),6)]./(PV(fluid).*U);
    end
    
    % Apparent viscosity and effective d_h based Re number for PL and YPPV
    Re_PL = Re_nom./eta_f_PL; Re_YPPV = Re_nom./eta_f_YPPV;
    % "Pipe" Reynolds number of Ouriemi et al. (2009, 2010) = superficial Re w/o solids bed
    if ii==1
        Re_pipe_PL = Re_PL;
        Re_pipe_YPPV = Re_YPPV;
    end
    Re_cr_PL = [3250-1150*n_P 3250-1150*n_A];
    % Bingham Reynolds numbner based on Bingham PV
    Re_Bi = Re_nom./my;
    He = Re_Bi.*Yi;
    Re_cr_Bi = CompCrReBi( He, length(PV) );
    
    % Loop all particle diameters
    % jj = 2; jj = 3;
    for jj = 1:length(d_p)
        
        %Particle Re number
        Re_p = rho_f*U.*d_p(jj)./eta_f_PL; Re_p(:,1:2,:)=0;

        % Critical Shields number (Soulsby 1997)
        d_star = (9.81.*(rho_s./rho_f-1)./(eta_f_PL./rho_f).^2).^(1/3).*d_p(jj); d_star(:,1:2,:)=0;
        Theta_cr = 0.24./d_star+0.055.*(1-exp(-0.02.*d_star)); Theta_cr(:,1:2,:)=0;
        
        % Sh number
        Sh_PL = eta_f_PL.*SR_PL./Sh_denom(jj); Sh_PL(:,1:2,:)=0;
        Sh_YPPV = eta_f_YPPV.*SR_YPPV./Sh_denom(jj); Sh_YPPV(:,1:2,:)=0;

        % Normalize Sh number w/ critical Shields number
        Sh_PL = Sh_PL./Theta_cr;
        Sh_YPPV = Sh_YPPV./Theta_cr;
        
        
        % Ar number
        Ar_PL = Ar_nom(jj)/eta_f_PL.^2; Ar_PL(:,1:2,:)=0;
        Ar_YPPV = Ar_nom(jj)/eta_f_YPPV.^2; Ar_YPPV(:,1:2,:)=0;
        % Ar_PL .* ((h_f./d_p(jj)).^2).*ones(size(Ar_PL));
        
        PlotOuriemiFlowPatternMap;
        PlotReShMap;
        counter=counter+1;
        
    end
end


figure(fig_Ouriemi);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


figure(fig_Re_Sh);

% Print to files
fig_name = num2str(get(gcf,'Number'));
set(gcf,'PaperPositionMode','auto') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'
