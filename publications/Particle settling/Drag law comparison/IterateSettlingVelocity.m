%   Iteratively compute settling velocity fur current case


% Estimate Stokes settling velocity based on H2O [m/s]
v_set = (rho_p-rho_f).*9.81.*(Cases.d_p(j,1)).^2./18./0.00102;
v_set_old = 0;

% Shear rate
SR = v_set./(0.5.*Cases.d_p(j,1));

% Accuracy
accuracy = 1e-8;


% Iteratively compute settling velocity
while abs(v_set - v_set_old) > accuracy % eps(v_set)
    
    % Update settling velocity
    v_set_old = v_set;
    
    switch draglaw
        case 'SchillerNaumann1935'
            
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = SchillerNaumann1935(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
            
            % Save results in table
            SN1935.SR(j,1) = SR;
            SN1935.eta(j,1) = eta;
            SN1935.Re_p(j,1) = Re_p;
            SN1935.Re_p2(j,1)= Re_p_Standard;
            SN1935.c_D(j,1) = c_D;
            SN1935.v_set(j,1) = v_set;
            
            SN1935.index(j,1) = j;
        
            
        case 'Acharya1976'
                
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = Acharya1976(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
            
            % Save results in table
            A1976.SR(j,1) = SR;
            A1976.eta(j,1) = eta;
            A1976.Re_p(j,1) = Re_p;
            A1976.Re_p2(j,1)= Re_p_Standard;
            A1976.c_D(j,1) = c_D;
            A1976.v_set(j,1) = v_set;
            
            A1976.index(j,1) = j;
        
            
       case 'ChhabraUhlherr1980'
                
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = ChhabraUhlherr1980(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Carreau.mu_inf{k},...
                Carreau.mu_0{k},...
                Carreau.lambda{k},...
                Carreau.n{k},...
                rho_f,...
                rho_p);
            
            % Save results in table
            CU1980.SR(j,1) = SR;
            CU1980.eta(j,1) = eta;
            CU1980.Re_p(j,1) = Re_p;
            CU1980.Re_p2(j,1)= Re_p_Standard;
            CU1980.c_D(j,1) = c_D;
            CU1980.v_set(j,1) = v_set;
            
            CU1980.index(j,1) = j;
        
            
        case 'Shah1986'
                       
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = Shah1986(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
                        
            % Save results in table
            S1986.SR(j,1) = SR;
            S1986.eta(j,1) = eta;
            S1986.Re_p(j,1) = Re_p;
            S1986.Re_p2(j,1)= Re_p_Standard;
            S1986.c_D(j,1) = c_D;
            S1986.v_set(j,1) = v_set;
            
            S1986.index(j,1) = j;
            
            
        case 'Ceylan1999'
                       
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = Ceylan1999(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
                        
            % Save results in table
            C1999.SR(j,1) = SR;
            C1999.eta(j,1) = eta;
            C1999.Re_p(j,1) = Re_p;
            C1999.Re_p2(j,1)= Re_p_Standard;
            C1999.c_D(j,1) = c_D;
            C1999.v_set(j,1) = v_set;
            
            C1999.index(j,1) = j;
            
            
            
        case 'Renaud2004'
                        
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = Renaud2004(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
                                   
            % Save results in table
            R2004.SR(j,1) = SR;
            R2004.eta(j,1) = eta;2
            R2004.Re_p(j,1) = Re_p;
            R2004.Re_p2(j,1)= Re_p_Standard;
            R2004.c_D(j,1) = c_D;
            R2004.v_set(j,1) = v_set;
            
            R2004.index(j,1) = j;
            
        case 'Shah2007'
                       
            % Compute settling velocity
            [SR,...
                eta,...
                Re_p,...
                Re_p_Standard,...
                c_D,...
                v_set] = Shah2007(v_set,...
                Cases.d_p(j,1),...
                SR,...
                Cross.mu_inf{k},...
                Cross.mu_0{k},...
                Cross.lambda{k},...
                Cross.n{k},...
                rho_f,...
                rho_p);
                        
            % Save results in table
            S2007.SR(j,1) = SR;
            S2007.eta(j,1) = eta;
            S2007.Re_p(j,1) = Re_p;
            S2007.Re_p2(j,1)= Re_p_Standard;
            S2007.c_D(j,1) = c_D;
            S2007.v_set(j,1) = v_set;
            
            S2007.index(j,1) = j;
            
        otherwise
            warning('Unexpected drag law. No plot created.');
    end
    
end % of while
            
            