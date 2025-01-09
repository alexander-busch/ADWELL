clc;
close all;
clear all;
addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange';

fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-12 - On the validity of the two-fluid-ktgf approach for dense gravity-driven granular flows\Figures\';

% Size and symbols for markers representing particle diameters    
markersizelist = [30,80,200];
markersymbolslist = {'o','d','s'};

% Colors for different fluids
markercolorlist = {'[0 0 1]' 'none'; ... % Air
    '[0 0 1]' '[0 0 1]'; ... H2O
    '[0 1 0]' '[0 1 0]'; ... PAC3
    '[1 0 0]' '[1 0 0]'}; % PAC4

% Load CFD results
load('results.mat');

% List of initial solid volume fractions investigated
alpha_s0_list = unique(AllCases.alpha_s_0); % [0.55 0.59 0.6]; % alpha_s0 = 0.55;

% List of fluids investigated
fluid_list = unique(AllCases.Fluid); % {'air' 'h2o' 'pac2' 'pac4'};

% List of aspect ratios investigated
AR_list = unique(AllCases.AR); % [2 3];

% Array of case_list investigated (combinations of cliff and particle scales)
case_list = unique(AllCases.Index); % 1:1:9;

% Loop initial solid volume fraction
for aa = 1:length(alpha_s0_list)
    % aa = 3
    
    % Loop aspect ratios
    for bb = 1:length(AR_list)
        % bb = 3
        
        % Loop fluids
        for cc = 1:length(fluid_list)
            % cc = 4

            % Get current variables and create "case2print" table
            alpha_s0 = alpha_s0_list(aa);
            AR = AR_list(bb);
            fluid = fluid_list{cc};
            case_list2print = AllCases((AllCases.alpha_s_0==alpha_s0) & (AllCases.AR==AR) & strcmp(AllCases.Fluid,fluid),:);

            % Create dimensional deposit height vs run-out distance figure
            % w/ 9 subplots
            fig_title = ['y vs x, \alpha_s_,_0 = ' num2str(alpha_s0) ', aspect ratio a = ' num2str(AR) ', ' fluid ];
            fig_handle = CreateFigure( fig_title, '', '', 'off', [1 1 29.7 15]);
            fig_sub = tight_subplot(3,3,[.1 .06],[.06 .03],[.06 .02]);

            % Reposition figure on screen
            set_fig_position( gcf, 1 );

            % Create non-dimensional deposit height vs aspect ratio and run-out
            % distance vs aspect ratio figures
            if bb==1            
                [fig_AR_x{cc}, fig_AR_y{cc}, fig_Shape_y_0{cc}] = CreateScalingFigures(fluid, alpha_s0);
                % close(fig_Shape_y_0{cc});
            end

            % Reposition figures on screen
            set_fig_position( fig_AR_x{cc}, 2 );
            set_fig_position( fig_AR_y{cc}, 3 );
            set_fig_position( fig_Shape_y_0{cc}, 4 );

            % Loop case_list (combinations of cliff and particle scales)
            for dd = 1:length(case_list)
                % bb = 2
                % bb = 3

                % Select dimensional deposit height vs run-out distance figure
                % and subplot for current case and plot dimensional data
                figure(fig_handle);
                subplot(fig_sub(dd));
                hold on;
                axis equal;
                box on;
                grid on;
                set(gca,...
                    'XScale','lin',...
                    'YScale','lin',...
                    'xlim', [0 case_list2print.x_max(dd)],...
                    'ylim', [0 case_list2print.y_max(dd)],...
                    'FontSize',10);

                % Plot IC as specified
                txt='k-';
                plot([0 case_list2print.x_0(dd)],[case_list2print.y_0(dd) case_list2print.y_0(dd)],txt);
                plot([case_list2print.x_0(dd) case_list2print.x_0(dd)],[0 case_list2print.y_0(dd)],txt);

                % Get all shapes for current case, loop them and plot them
                shapes_x = case_list2print.shapes_x{dd,:};
                shapes_y = case_list2print.shapes_y{dd,:};
                for ee = 1:min(size(shapes_x))
% figure
                    % Plot first shape with dotted lines as it is the actual IC
                    if (ee==1 && case_list2print.alpha_s_0(dd)==0.59)
                        plot([0 case_list2print.x_0(dd)],[case_list2print.y0_true(dd) case_list2print.y0_true(dd)],'k--');
                        plot([case_list2print.x_0(dd) case_list2print.x_0(dd)],[0 case_list2print.y0_true(dd)],'k--');
                    else
                        plot(shapes_x(:,ee),shapes_y(:,ee),'-k');
                    end

                end % of loop shapes

                % Select marker type representing scale of cliff
                if dd<4
                    markersymbol = markersymbolslist{1};
                elseif dd<7
                    markersymbol = markersymbolslist{2}; 
                else
                    markersymbol = markersymbolslist{3};
                end

                % Select marker size representing scale of particle diameter
                if case_list2print.d_s(dd)==case_list2print.d_s(1)
                    markersize = markersizelist(1);
                elseif case_list2print.d_s(dd)==case_list2print.d_s(2)
                    markersize = markersizelist(2);
                else
                    markersize = markersizelist(3);
                end

                % Mark the final run-out length and the final deposit height
                scatter(case_list2print.x_f_CFD(dd),0,markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');
                scatter(min(xlim),case_list2print.y_f_CFD(dd),markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');              
                
                % Compute Froude number
                switch case_list2print.Fluid{dd}
                    case 'air'
                        mu_0 = 1.79e-5;
                    case 'h2o'
                        mu_0 = 1.002e-3;
                    case 'pac2'
                        % Simply use Newtionian low-Shear limit
                        mu_0 = 7.21e-2;
                        mu_inf = 1.002e-3;
                        lambda_Cr = 1.090e-2;
                        n_Cr = 0.586;
                    case 'pac4'
                        mu_0 = 2.14e-1; 
                        mu_inf = 1.002e-3;
                        lambda_Cr = 2.610e-2;
                        n_Cr = 0.608;
                    otherwise
                       mu_0 = 2.14e-1;
                       warning('Unexpected eta. No Fr created.');
                end                   
                
                % Estimate Stokes settling velocity based on zero-shear-viscosity  [m/s]
                v_set = (case_list2print.rho_s(dd)-case_list2print.rho_f(dd))*9.81*(case_list2print.d_s(dd))^2/18/mu_0;
                
                % If Non Newtonian compute settling velocity iteratively
                if any( ismember(case_list2print.Fluid{dd},{'pac2', 'pac4'}))               

                    v_set_old = 0;

                    % Accuracy
                    accuracy = 1e-8;          

                    % Iteratively compute settling velocity
                    while abs(v_set - v_set_old) > accuracy % eps(v_set)

                        % Update settling velocity
                        v_set_old = v_set;

                        % Shear rate
                        SR = v_set/case_list2print.d_s(dd);

                        % Apparent viscosity
                        eta_f = mu_inf+(mu_0-mu_inf)/(1+(lambda_Cr*SR)^(1-n_Cr));

                        % Particle Reynolds number 
                        Re_p = case_list2print.rho_f(dd)*v_set*case_list2print.d_s(dd)/eta_f;

                        % Schiller & Naumann (1935) drag coefficient
                        c_D = (24/Re_p)*(1+0.15*(Re_p^0.687));

                        % Schiller & Naumann (1935) high Re correction
                        c_D(Re_p > 1000)=0.44;

                        % Settling velocity
                        v_set = (4*case_list2print.d_s(dd)*9.81/3/c_D*(case_list2print.rho_s(dd)/case_list2print.rho_f(dd)-1))^0.5;

                    end
                end

                Fr = case_list2print.rho_s(dd)*(case_list2print.d_s(dd))^3*9.81*(case_list2print.rho_s(dd)-case_list2print.rho_f(dd))/(18^2*mu_0);
                I_square = v_set^2/9.81/case_list2print.y_0(dd);
                
                % Add "simulation ongoing" text or case 1-9 text
                if isnan(case_list2print.y_f_CFD(dd))
                    %text(0.5*max(xlim),0.5*max(ylim),['Simulation ongoing (' date ')'],'fontsize',14,'color','r','HorizontalAlignment','center');
                else
                    % text(0.5*max(xlim),0.5*max(ylim),['Case ' num2str(dd) ', Fr = ' num2str(Fr,2) ', I^2 = ' num2str(I_square,2)],'fontsize',10,'color','k','HorizontalAlignment','center');
                    text(0.5*max(xlim),0.5*max(ylim),['Case ' num2str(dd)],'fontsize',14,'color','k','HorizontalAlignment','center');
                    % {['Case ' num2str(dd)];['d_p = ' num2str(case_list2print.d_s(dd)*1000) ' mm, x_0 = ' num2str(case_list2print.x_0(dd)) ' m' ]}
                end
                
                
                
                % Plot on non-dimensional shape figure
                % Scaling add. to Lube et al. (2005)
                figure(fig_Shape_y_0{cc});
                x_n = shapes_x(:,end)/case_list2print.x_f_CFD(dd);
                y_n = shapes_y(:,end)/case_list2print.y_f_CFD(dd);
                scatter(x_n,smooth(y_n),'.'); % Plot IF

                % Normalize dimensional run-out length and deposit height
                x_n = (case_list2print.x_f_CFD(dd)-case_list2print.x_0(dd))/case_list2print.x_0(dd);
                % x_n = case_list2print.x_f_CFD(dd)/case_list2print.x_0(dd);
                % y_n = case_list2print.y_f_CFD(dd)/case_list2print.x_0(dd); % Normalized with x0 add. to Lube et al. (2005)
                 y_n = case_list2print.y_f_CFD(dd)/case_list2print.y_0(dd); % Normalized with y0
               % y_n = case_list2print.y_f_CFD(dd)/case_list2print.y0_true(dd); % Normalized with true y0 (result of pre-simulation

                % Plot on non-dimensional f(AR) figures
                figure(fig_AR_x{cc});
                 scatter(case_list2print.AR(dd),x_n,markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');
                %scatter(case_list2print.AR_true(dd),x_n,markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');
                figure(fig_AR_y{cc});
                scatter(case_list2print.AR(dd),y_n,markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');
                % scatter(case_list2print.AR_true(dd),y_n,markersize,markersymbol,'MarkerEdgeColor',markercolorlist{cc},'MarkerFaceColor','none');

            end % of loop case_list
                            
            % Export dimensional shape figure -----------------------------
            figure(fig_handle);
            
            % Set paper position mode for printing
            set(gcf,'PaperPositionMode','auto');
            fig_handle_pos = fig_handle.PaperPosition;
            fig_handle.PaperSize = [fig_handle_pos(3) fig_handle_pos(4)];
            fig_handle.PaperUnits = 'centimeters';

            % Save Figure to File Format
            fig_name = ['alpha_s0=' num2str(alpha_s0) '_a=' num2str(AR) '_' fluid '_shape_dimensional'];
            print(gcf,[fig_path fig_name '.pdf'],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
            print(gcf,[fig_path fig_name '.svg'],'-dsvg'); 
            
            
            % Export non-dimensional shape figure -------------------------
%             figure(fig_Shape_y_0{cc});
% 
%             % Set paper position mode for printing
%             set(gcf,'PaperPositionMode','auto');
%             fig_handle_pos = fig_handle.PaperPosition;
%             fig_handle.PaperSize = [fig_handle_pos(3) fig_handle_pos(4)];
%             fig_handle.PaperUnits = 'centimeters';
%             
%             % Save Figure to File Format
%             fig_name = ['alpha_s0=' num2str(alpha_s0) '_a=' num2str(AR) '_' fluid '_shape_non-dimensional'];
%             print(gcf,[fig_path fig_name '.pdf'],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
%             print(gcf,[fig_path fig_name '.svg'],'-dsvg'); 
            
            % Export AR figure --------------------------------------------
            if AR==AR_list(end)

                figure(fig_AR_x{cc});

                % Set paper position mode for printing
                set(gcf,'PaperPositionMode','auto');
                fig_handle_pos = fig_AR_x{cc}.PaperPosition;
                fig_AR_x{cc}.PaperSize = [fig_handle_pos(3) fig_handle_pos(4)];
                fig_AR_x{cc}.PaperUnits = 'centimeters';

                % Save Figure to File Format
                fig_name = ['alpha_s0=' num2str(alpha_s0) '_' fluid '_scaling_x'];
                print(gcf,[fig_path fig_name '.pdf'],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
                print(gcf,[fig_path fig_name '.svg'],'-dsvg'); 


                figure(fig_AR_y{cc});

                % Set paper position mode for printing
                set(gcf,'PaperPositionMode','auto');
                fig_handle_pos = fig_AR_y{cc}.PaperPosition;
                fig_AR_y{cc}.PaperSize = [fig_handle_pos(3) fig_handle_pos(4)];
                fig_AR_y{cc}.PaperUnits = 'centimeters';
                
                % Save Figure to File Format
                fig_name = ['alpha_s0=' num2str(alpha_s0) '_' fluid '_scaling_y'];
                print(gcf,[fig_path fig_name '.pdf'],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
                print(gcf,[fig_path fig_name '.svg'],'-dsvg'); 
            end

        end % of loop fluids
        
    end % of loop aspect ratio

end % of loop initial solid volume fraction