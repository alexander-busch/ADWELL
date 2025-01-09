function [ Re_p, c_D, v_set ] = IterateSettlingVelocity( d_p, rho_p, rho_f, mu_inf, mu_0, lambda_Cr, n_Cr, SR_mean, SR_fluc )
%   IterateSettlingVelocity Iteratively compute settling velocity
%   Detailed explanation goes here

	
% Estimate Stokes settling velocity based on H2O [m/s]
v_set = (rho_p-rho_f).*9.81.*(d_p).^2./18./0.00102;

% Initialize v_set_old
v_set_old = zeros(length(v_set),1);

% Initialize eta
eta = 0.00102;

% Iteratively compute settling velocity
while (max(abs(v_set-v_set_old)) >= 1e-4)
    
    % Update settling velocity
    v_set_old = v_set;


    

    % Schiller-Naumann (1933)----------------------------------------------

    SR_p = v_set./(0.5.*d_p);
    
    SR = (SR_fluc.^2 + SR_mean.^2 + SR_p.^2).^0.5;
    
    eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);
    
    Re_p = Rep_Standard(rho_f, v_set, eta, d_p);
        
    c_D = cD_SchillerNaumann1935(Re_p);
    
    % ---------------------------------------------------------------------
    
    
    
    % Variant 2 -----------------------------------------------------------
    % 1) Effective Cross Shear Rate --> Effective viscosity
    % 2) Particle Reynolds number based on effective viscosity
    % 3) Schiller-Naumann drag law

%     c1 = -((1-n_Cr).*(mu_0.^2).*SR.^(1-n_Cr))./(lambda_Cr.*(1+lambda_Cr.*SR.^(1-n_Cr)).^2);
%     c2 = mu_0./(1+lambda_Cr.*SR.^(1-n_Cr));
%     c3 = mu_0.*SR./(1+lambda_Cr.*SR.^(1-n_Cr));
%     n_dash = SR.*(c1+c2)./c3;
%     SR = (3.*n_dash+1)./(4.*n_dash).*SR;
%     
%     eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);
%     
%     Re_p = Rep_Standard(rho_f, v_set, eta, d_p);
%     
%     c_D = cD_SchillerNaumann1935(Re_p);
    
    % Gives almost identical results as Variant 1 but cannot handle case of
    % water with lambda_Cr = 0
    
    % End of Variant 2 ----------------------------------------------------
     
    
    
    % Renaud (2004) -------------------------------------------------------

%     [ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );
%         
%     SR_p = SR_Renaud2004( v_set, n_PL, d_p )
%  
%     SR = (SR_fluc.^2 + SR_mean.^2 + SR_p.^2).^0.5;
%     
%     Re_p = Rep_Renaud2004(rho_f, v_set, n_PL, K_PL, d_p);
%     
%     c_D = cD_Renaud2004(Re_p, n_PL);

    % --------------------------------------------------------------------- 

    
    
    % Shah (2007) ---------------------------------------------------------
    
%     SR_p = v_set./(0.5.*d_p);
% 
%     SR = (SR_fluc.^2 + SR_mean.^2 + SR_p.^2).^0.5;
%     
%     [ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR );
%     
%     Re_p = Rep_Shah2007(rho_f, v_set, n_PL, K_PL, d_p);
%     
%     c_D = cD_Shah2007(Re_p, n_PL);

    % --------------------------------------------------------------------- 
    
    
    
    % Shah & Uhlherr (1980)------------------------------------------------

%     Re_p = rho_f.*v_set.*d_p./mu_0;
%     Ca = 2*lambda_Ca.*v_set./d_p;
%     
%     c_D = 24./Re_p.*(1 + 0.15.*Re_p.^0.687).*(1+0.65.*(n_Ca-1)*Ca.^0.2);
    
    % --------------------------------------------------------------------- 
    
    
    
    
    
    % Settling velocity
    v_set = (4.*d_p.*9.81./3./c_D.*(rho_p./rho_f-1)).^0.5;
end % of while

% Final Particle Reynolds number
Re_p = Rep_Standard(rho_f, v_set, eta, d_p);

% Final Drag Coefficient
c_D = 4.*d_p.*9.81.*(rho_p./rho_f-1)./3./v_set.^2;
c_D = cD_SchillerNaumann1935(Re_p);

end % of function

