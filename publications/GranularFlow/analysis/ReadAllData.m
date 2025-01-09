clc;
close all;
clear all;
addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange';

filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive\';


% List of initial solid volume fractions investigated
alpha_s0list = [0.55 0.59 0.6]; % alpha_s0 = 0.55; % To be changed later

% List of aspect ratios investigated
ARlist = [1 2 3 ];

% List of fluids investigated
fluidlist = {'air' 'h2o' 'pac2' 'pac4'};

% Read metadata
% [num,txt,raw] = xlsread([filepath 'metadata.xlsx'],'metadata','B2:D100');
% metadata = readtable([filepath 'metadata.xlsx'],'Sheet','metadata','Range','B1:D100', 'ReadVariableNames', true);

% Create intermediate case table
index = 1:1:9;
Cases = table;
Cases.Index = index';                           % Table index
Cases.x_max = repelem([0.10 1.0 10]',3);        % x dimension of computational domain
Cases.y_max = repelem([0.04 0.4 04]',3);        % y dimension of computational domain
Cases.dx = repelem([0.002 0.02 0.2]',3);        % Mesh size
Cases.d_s = repmat([0.0001 0.001 0.01]',3,1);	% Solid particle diameter
Cases.rho_s = repmat(2650,length(index),1);     % Solid density
Cases.theta_s = repmat(45,length(index),1);     % Solid angle of internal friction
Cases.rho_f = repmat(1.225,length(index),1);    % Fluid density
Cases.x_0 = 0.1*Cases.x_max;                    % Initial width
Cases.y_0 = Cases.y_max-Cases.x_0;              % Initial height
% Cases.alpha_s_0 = 0.61*ones(length(index),1);   % Initial average volume fraction


% Loop initial solid volume fraction
for aa = 1:length(alpha_s0list)
    % aa = 3
    
    % Set initial average volume fraction
    Cases.alpha_s_0 = alpha_s0list(aa)*ones(length(index),1);
    
    % Loop aspect ratio
    for bb = 1:length(ARlist)
        % bb = 3
        
        % Loop fluids
        for cc = 1:length(fluidlist)
            % cc = 4
            
            % Get current variables and write to intermediate "case" table
            AR = ARlist(bb);
            Cases.AR(:,1)=AR;
            
            fluid = cell(9,1); [fluid{:}] = deal(fluidlist{cc});
            Cases.Fluid(:,1) = fluid;
            fluid = fluidlist{cc}; 
            
            Cases.y_0(:,1)=Cases.x_0(:,1).*AR;
            
            % Set number of CFD-Post export data files per case and fluid density
            if alpha_s0list(aa)==0.59
                % Cases with presettling simulations
                if strcmp(fluid,'air')
                    cfdpostexport = 4; % Number of exported data files per case, air 1-4, liquids 1-10
                    Cases.rho_f = repmat(1.225,length(index),1);    % Fluid density
                else
                    cfdpostexport = 10; % Number of exported data files per case, air 1-4, liquids 1-10
                    if strcmp(fluid,'h2o')
                        Cases.rho_f = repmat(998,length(index),1);    % Fluid density
                    else
                        Cases.rho_f = repmat(1000,length(index),1);    % Fluid density
                    end
                end
            elseif alpha_s0list(aa)==0.6
                % Cases with patch, alpha_S = 0.6
                if strcmp(fluid,'air')
                    cfdpostexport = 2; % Number of exported data files per case
                    Cases.rho_f = repmat(1.225,length(index),1);    % Fluid density
                elseif strcmp(fluid,'h2o')
                    cfdpostexport = 6; % Number of exported data files per case
                    Cases.rho_f = repmat(998,length(index),1);    % Fluid density
                elseif strcmp(fluid,'pac2')
                    cfdpostexport = 12; % Number of exported data files per case
                    Cases.rho_f = repmat(998,length(index),1);    % Fluid density
                else
                     cfdpostexport = 13; % Number of exported data files per case
                    Cases.rho_f = repmat(1000,length(index),1);    % Fluid density
                end
            else
                % Cases with patch, alpha_S = 0.55
                if strcmp(fluid,'air')
                    cfdpostexport = 2; % Number of exported data files per case
                    Cases.rho_f = repmat(1.225,length(index),1);    % Fluid density
                elseif strcmp(fluid,'h2o')
                    cfdpostexport = 2; % Number of exported data files per case
                    Cases.rho_f = repmat(998,length(index),1);    % Fluid density
                elseif strcmp(fluid,'pac2')
                    cfdpostexport = 4; % Number of exported data files per case
                    Cases.rho_f = repmat(998,length(index),1);    % Fluid density
                else
                     cfdpostexport = 6; % Number of exported data files per case
                    Cases.rho_f = repmat(1000,length(index),1);    % Fluid density
                end
            end

            % Loop all cases
            for dd = 1:length(index)
                % Debugging
                % bb = 2
                % bb = 3

                % Loop all CFD-Post export files for current case
                for ee = 0:cfdpostexport
                    % ee = 1
                    % ee = 2
                    % ee = 3
                    % ee = 4
                    % ee = 10

                    
%                     if AR==2
%                         [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports\a=' num2str(AR) '\' fluid '\case' num2str(Cases.Index(dd)) '_t' num2str(ee,'%02d') '.csv']);     
%                     else
% 
%                         if cfdpostexport<10
%                             [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports\a=' num2str(AR) '\' fluid '\case' num2str(Cases.Index(dd)) '_t' num2str(ee) '.csv']);     
%                         else
%                             [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports\a=' num2str(AR) '\' fluid '\case' num2str(Cases.Index(dd)) '_t' num2str(ee,'%02d') '.csv']);     
%                         end
% 
%                     end

                    % Determine numebr format of CFD post export files
                    if (alpha_s0list(aa)==0.59 && AR==3 && strcmp(fluid,'air'))
                        cfdpostexport_num = num2str(ee);
                    else
                        cfdpostexport_num = num2str(ee,'%02d');
                    end
                    
                    % Read cfd post data
                    cfdpostexport_filename = [pwd '\cfdpostexports\alpha_s=' num2str(alpha_s0list(aa)) '\a=' num2str(AR) '\' fluid '\case' num2str(Cases.Index(dd)) '_t' cfdpostexport_num '.csv'];
                    if exist(cfdpostexport_filename)
                        [x,y,vof] = importCFDpostdata(cfdpostexport_filename);

                        % surf(x,y,vof)
                        % Shift x-data such that x0 = 0
                        x = x + Cases.x_max(dd)/2;

                        % Transpose vof and remove zeros
                        vof=vof';
                        vof(vof==0)=nan;

                        % Determine gradient
                        [dvofdx,dvofdy] = gradient(vof);
                        grad_mag = (dvofdx.^2+dvofdy.^2).^0.5;

                        % Determine interface (IF) y = f(x)
                        IF = zeros(length(x),3);

                        % Loop all x and evaluate y-direction vectors 
                        for ff = 1:length(x)
                            %ee=65
                            IF(ff,1) = x(ff); %+Cases.x_max(bb)/2; % Create IF profile: x

                            % Current vof vector in y-direction
                            current_vof = vof(:,ff);

                            if ((max(current_vof)<0.5) || (isnan(max(current_vof))))
                                IF(ff,2) = nan;
                            else
                                IF(ff,2) = max(y(current_vof>0.5));
                            end

                            % Current grad_mag vector in y-direction
                            current_grad_mag = grad_mag(:,ff);
                            if isnan(max(current_grad_mag))
                                IF(ff,3) = nan;
                            else
                                IF(ff,3) = y(current_grad_mag==max(current_grad_mag)); % Create IF profile: y      
                            end
                        end % of loop all x and evaluate y-direction vectors

                        % Indices of unique values in IF(:,1) = x    
                        [~, ind] = unique(IF(:,1), 'rows');
                        % Interpolate IF in order to close the data gap @ x = 0.05, 0.05, 5
                        xi = linspace(0,Cases.x_max(dd));
                        yi = interp1(IF(ind,1),IF(ind,3),xi);

                        % Number of cells for width of initial conditions case
                        cells = round(Cases.x_0(dd)/Cases.dx(dd));

                        % Final deposit height
                        Cases.y_f_CFD(dd) = max(yi(1:cells));

                        % Remove wrong data points at end
                        logic = yi>Cases.y_f_CFD(dd);
                        yi(logic) = nan;

                        % Final run-out
                        if dd~=3
                            logic = yi>=Cases.d_s(dd); % Get logical array which contains zeros for all positions where the deposit height is less than the particle diameter
                            % Check if logical array contains zeros
                            if nnz(logic)==length(logic)
                                Cases.x_f_CFD(dd) = xi(end); 
                            else
                                position = find(logic==0, 1, 'first');           
                                if logic(position+1) ~= 0
                                    position = find(logic==0, 1, 'last');
                                end
                                Cases.x_f_CFD(dd) = xi(position);
                            end
                        else
                            position = find(isnan(yi), 1, 'first')-1;
                            Cases.x_f_CFD(dd) = xi(position);
                        end

                        % Save shapes to matrix
                        if ee==0
                            shapes_x = xi';
                            shapes_y = yi';
                        else
                            shapes_x = [shapes_x xi'];
                            shapes_y = [shapes_y yi'];
                        end

                    else
                        % No data for current case
                        Cases.y_f_CFD(dd) = nan;
                        Cases.x_f_CFD(dd) = nan;
                        shapes_x = nan;
                        shapes_y = nan;
                    end
                end % of loop cfdpostexport

                Cases.shapes_x{dd} = shapes_x;
                Cases.shapes_y{dd} = shapes_y;

                % Determine true y0 and AR
                Cases.y0_true(dd) = nanmean(shapes_y(:,1));
                Cases.AR_true(dd) = Cases.y0_true(dd)/Cases.x_0(dd);

            end % of loop cases

            % Assemble final table       
            if ((aa==1) && (bb==1)&&(cc==1))
                AllCases = Cases;
            else
                AllCases = [AllCases; Cases];
            end
            
        end % of loop fluids
        
    end % of loop aspect ratio
    
end % of loop initial solid volume fraction



save('results.mat','AllCases');