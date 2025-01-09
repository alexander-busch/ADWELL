% Read exp results and write time, x, y, Re_p vectors to table


if i == 1 % H2O
    if d_p(j) == 0.00116
        if U(k) == 0.048
%             EXP.t(row,1) = {0};
%             EXP.X(row,1) = {0};
%             EXP.Y(row,1) = {0};
%             EXP.vx0(row,1) = (0);
%             EXP.vy0(row,1) = (0);
        elseif U(k) == 0.085
            
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'A2:K150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m
            
            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
                EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.16;
            end
            
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison Water and 3DCFD','A5:G150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            CFD_Khatibi.vx0(row,1) = data(1,6);
            CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.16;
            
        else
        end
    elseif d_p(j) == 0.002
        if U(k) == 0.048
%             EXP.t(row,1) = {0};
%             EXP.X(row,1) = {0};
%             EXP.Y(row,1) = {0};
%             EXP.vx0(row,1) = (0);
%             EXP.vy0(row,1) = (0);
        elseif U(k) == 0.085
                        
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'R2:AB150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 

            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
            	EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.28;
            end
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison Water and 3DCFD','J5:P150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            CFD_Khatibi.vx0(row,1) = data(1,6);
            CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.28;
        else
        end
    else
    end

elseif i == 2 % PAC2
    if d_p(j) == 0.002
        if U(k) == 0.048
                                    
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'A2:K150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 

            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
            	EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.02;
            end                                    
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 2 and 3DCFD','J5:N150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.02;
            
        elseif U(k) == 0.085
%             EXP.t(row,1) = {0};
%             EXP.X(row,1) = {0};
%             EXP.Y(row,1) = {0};
%             EXP.vx0(row,1) = (0);
%             EXP.vy0(row,1) = (0);
        else
        end
    elseif d_p(j) == 0.003
        if U(k) == 0.048
                                                
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'R2:AB150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 
            
            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
            	EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.05;
            end                                                
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 2 and 3DCFD','X5:AA150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.05;
            
        elseif U(k) == 0.085
%             EXP.t(row,1) = {0};
%             EXP.X(row,1) = {0};
%             EXP.Y(row,1) = {0};
%             EXP.vx0(row,1) = (0);
%             EXP.vy0(row,1) = (0);
        else
        end
    else
    end
    
    
    
elseif i==3 % PAC4

    if d_p(j) == 0.002
        if U(k) == 0.048
                                                            
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'AI2:AS150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 

            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
            	EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.007;
            end
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 4 and 3DCFD','AO4:AR150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.007;
            
            
        elseif U(k) == 0.085
                                                            
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'A2:K150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 
            
            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
            	EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.007;
            end
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 4 and 3DCFD','I5:L150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.007;
            
        else
        end
    elseif d_p(j) == 0.003
        if U(k) == 0.048
                                                            
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'AX2:BH150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 
            
            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
                EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.015;
            end                                    
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 4 and 3DCFD','BD4:BG150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.015;
           
        elseif U(k) == 0.085
                                                            
            % Khatibi et al. (2016) - Experimental data
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', fluidlist{i},'Q2:AA150');
            EXP.t(row,1) = {data(:,1)};
            EXP.X(row,1) = {data(:,3)/1000}; % mm --> m 
            EXP.Y(row,1) = {data(:,5)/1000}; % mm --> m 
            
            if v_ini == 0
                EXP.vx0(row,1) = data(1,9);
                EXP.vy0(row,1) = data(1,11);
            else
                EXP.vx0(row,1) = 0.005;
                EXP.vy0(row,1) = -0.015;
            end
                                                
            
            % Khatibi et al. (2016) - CFD data            
            data = xlsread('exp\x and y particle trajectory - AS and AB.xlsx', 'Comparison PAC 4 and 3DCFD','Y5:ABG150');
            CFD_Khatibi.t(row,1) = {data(:,1)};
            CFD_Khatibi.X(row,1) = {data(:,3)/1000}; % mm --> m 
            CFD_Khatibi.Y(row,1) = {data(:,4)};
            % CFD_Khatibi.vx0(row,1) = data(1,6);
            % CFD_Khatibi.vy0(row,1) = data(1,7);
            
            CFD_Khatibi.vx0(row,1) = 0.005;
            CFD_Khatibi.vy0(row,1) = -0.015;
           
        else
        end
    else
    end
else
end
