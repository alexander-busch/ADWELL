clc;
% close all;
clear all;
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';

% Size for symbols representing particle diameters    
markersizelist = [30,80,200];

CreateScalingFigures;

Cases = table;
Cases.x_max = [0.10 1.0 10]';       % x dimension of computational domain
Cases.y_max = [0.04 0.4 04]';       % y dimension of computational domain
Cases.dx = [0.002 0.02 0.2]';       % Mesh size
Cases.d_s = [0.0001 0.001 0.01]';   % Solid particle diameter
Cases.rho_s = 2650*ones(3,1);       % Solid density
Cases.theta_s = 45*ones(3,1);       % Solid angle of internal friction
Cases.rho_f = 1.225*ones(3,1);      % Fluid density
Cases.x_0 = [0.01 0.1 01]';         % Initial width
Cases.y_0 = [0.03 0.3 03]';         % Initial height
Cases.alpha_s_0 = [0.5 0.5 0.5]';   % Initial average volume fraction
Cases = [Cases; Cases];
Cases.Name = {'small',...
    'medium',...
    'large',...
    'small_no_ps',...
    'medium_no_ps',...
    'large_no_ps'}';

for aa = 1:length(Cases.Name)
    % Debugging
    % aa = 2

    % Create dimensional deposit shape figure or select correct figure
    if aa<4
        switch aa
            case 1
                fig_small  = figure('color','w','Name','Dimensional deposit height as function of run-out distance (Shape of final profiles) - Small','Units','centimeters','Position',[10 10 12 8.45]);
                hold on;
                plot([0 0.01],[0.03 0.03],'k--');
                plot([0.01 0.01],[0 0.03],'k--');
            case 2
                fig_medium  = figure('color','w','Name','Dimensional deposit height as function of run-out distance (Shape of final profiles) - Medium','Units','centimeters','Position',[10 10 12 8.45]);
                hold on;
                plot([0 0.1],[0.3 0.3],'k--');
                plot([0.1 0.1],[0 0.3],'k--');
            case 3
                fig_large  = figure('color','w','Name','Dimensional deposit height as function of run-out distance (Shape of final profiles) - Large','Units','centimeters','Position',[10 10 12 8.45]);
                hold on;
                plot([0 1],[3 3],'k--');
                plot([1 1],[0 3],'k--');
            otherwise  
        end
        axis equal;
        box on;
        grid on;
        set(gca,...
            'XScale','lin',...
            'YScale','lin',...
            'xlim', [0 Cases.x_max(aa)],...
            'ylim', [0 Cases.y_max(aa)]);
    else
        switch aa
            case 4
                figure(fig_small);
            case 5
                figure(fig_medium);
            case 6
                figure(fig_large);
            otherwise  
        end
    end
    
    % Read cfd post data for t=0
    [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports_old\' Cases.Name{aa} '_0.csv']);
    % Number of cells for width of initial conditions case
    cells = round(Cases.x_0(aa)/Cases.dx(aa));
    % Iniialize looping variables
    index=1; IF = zeros(3,1);
    % Loop all x and y
    for ii = 1:cells
        for jj = 1:length(y)
            if vof(ii,jj)>0.1
                if vof(ii,jj)<0.55
                    % max(max(vof))
                    IF(1,index) = x(ii)+Cases.x_max(aa)/2; % Create IF profile: x
                    IF(2,index) = y(jj); % Create IF profile: y
                    IF(3,index) = vof(ii,jj);
                    index=index+1;
                else
                   
                end
            end
        end
    end  
    % x at t=0
    Cases.x_0(aa) = 0.1*Cases.x_max(aa);
    % y at t=0 
    Cases.y_0(aa) = mean(IF(2,1:cells))+ Cases.dx(aa)*mean(IF(3,1:cells));
    % Plot IC
    if aa<4
        txt='k-';
    else
        txt='b-';
    end
    plot([0 Cases.x_0(aa)],[Cases.y_0(aa) Cases.y_0(aa)],txt); % Plot IF
    plot([Cases.x_0(aa) Cases.x_0(aa)],[0 Cases.y_0(aa)],txt); % Plot IF
    
    % Effective aspect ratio for settled column
    Cases.AR(aa) = Cases.y_0(aa)/Cases.x_0(aa);    
    
    if aa==1
        for bb=1:4
            if bb==3
                
            else
                % Read cfd post data for final cliff-collapse
                [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports_old\' Cases.Name{aa} '_' num2str(bb) '.csv']);
                vof=vof';

                % Determine gradient
                [dvofdx,dvofdy] = gradient(vof);
                grad_mag = (dvofdx.^2+dvofdy.^2).^0.5;
            %     figure
            %     contour(x,y,grad_mag);
            %     contour(x,y,vof);
            %     hold on
            %     quiver(x,y,dvofdx,dvofdy)
            %     hold off

                % Set gradient to zero at domain boundaries
            %     for ii = 1:length(x)
            %         for jj = 1:length(y)
            %             if ((x(ii) <= -Cases.x_max(aa)/2+Cases.dx(aa)) || (y(jj) == 0))
            %                 grad_mag(jj,ii)=0;
            %             end
            %         end
            %     end

                for ii = 1:length(x)
                    for jj = 1:length(y)
                        if ((x(ii) <= -Cases.x_max(aa)/2) || (y(jj) == 0))
                            grad_mag(jj,ii)=0;
                        end
                    end
                end

            %     for ii = 1:length(x)
            %         for jj = 1:length(y)
            %             if grad_mag(jj,ii)<0.39
            %                 grad_mag(jj,ii)=0;
            %             end
            %         end
            %     end


                % Determine interface (IF) y = f(x)
                IF = zeros(3,1);
                index=2;
                % Loop all x and y
                for ii = 2:length(x)
                    for jj = 2:length(y)
                        if grad_mag(jj,ii)>0.9*max(max(grad_mag))
                            IF(1,index) = x(ii)+Cases.x_max(aa)/2; % Create IF profile: x
                            IF(2,index) = y(jj); % Create IF profile: y
                            index=index+1;
                        end
                    end
                end

                % Final deposit height
                IF(2,1)= max(IF(2,:));
                Cases.y_f_CFD(aa) = IF(2,1);

                % Final Run-out
                x_f = IF(2,:)>=Cases.d_s(aa); % Get logical array which contains zeros for all positions where the deposit height is less than the particle diameter
                % Check if logical array contains zeros
                if nnz(x_f)==length(x_f)
                    Cases.x_f_CFD(aa) = IF(1,end); 
                else
                    position = find(x_f==0, 1, 'first');
                    if x_f(position+1) ~= 0
                        position = find(x_f==0, 1, 'last');
                    end
                    Cases.x_f_CFD(aa) = IF(1,position);
                end

                % Plot IF on dimensional shape figure


                % Plot on dimensional shape figure
                switch aa
                    case 1
                        figure(fig_small);
                        markersize = markersizelist(1);
                        markertype = 'o';
                    case 2
                        figure(fig_medium);
                        markersize = markersizelist(2);
                        markertype = 'o';
                    case 3
                        figure(fig_large);
                        markersize = markersizelist(3);
                        markertype = 'o';
                    case 4
                        figure(fig_small);
                        markersize = markersizelist(1);
                        markertype = 'x';
                    case 5
                        figure(fig_medium);
                        markersize = markersizelist(2);
                        markertype = 'x';
                    case 6
                        figure(fig_large);
                        markersize = markersizelist(3);
                        markertype = 'x';
                    otherwise  
                end
                plot(IF(1,:),IF(2,:),'-k');
                scatter(Cases.x_f_CFD(aa),0,markertype,'b');
                scatter(min(xlim),Cases.y_f_CFD(aa),markertype,'g');

                % Plot on non-dimensional f(AR) figure
                figure(fig_AR);
                x_n = (Cases.x_f_CFD(aa)-Cases.x_0(aa))/Cases.x_0(aa);
                x_n = Cases.x_f_CFD(aa)/Cases.x_0(aa);
                yyaxis left;
                scatter(Cases.AR(aa),x_n,markersize,markertype,'b');
                y_n = Cases.y_f_CFD(aa)/Cases.x_0(aa); % Normalized with x0 acc. to Lube et al. (2005)
                y_n = y_n/Cases.AR(aa); % Normalized with y0
                yyaxis right;
                scatter(Cases.AR(aa),y_n,markersize,markertype,'g');

                % Plot on non-dimensional shape figure
                % Scaling acc. to Lube et al. (2005)
                figure(fig_Shape_y_0);
                x_n = IF(1,:)/Cases.x_f_CFD(aa);
                y_n = IF(2,:)/Cases.y_f_CFD(aa);
                if aa<4
                    scatter(x_n,smooth(y_n),markersize,'k'); % Plot IF
                else
                    % scatter(x_n,smooth(y_n),markersize,'b'); % Plot IF
                end
                
                
                
                
                
            end
            
        end
        
    else
        % Read cfd post data for final cliff-collapse
        [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports_old\' Cases.Name{aa} '_4.csv']);
        vof=vof';

        % Determine gradient
        [dvofdx,dvofdy] = gradient(vof);
        grad_mag = (dvofdx.^2+dvofdy.^2).^0.5;
    %     figure
    %     contour(x,y,grad_mag);
    %     contour(x,y,vof);
    %     hold on
    %     quiver(x,y,dvofdx,dvofdy)
    %     hold off

        % Set gradient to zero at domain boundaries
    %     for ii = 1:length(x)
    %         for jj = 1:length(y)
    %             if ((x(ii) <= -Cases.x_max(aa)/2+Cases.dx(aa)) || (y(jj) == 0))
    %                 grad_mag(jj,ii)=0;
    %             end
    %         end
    %     end

        for ii = 1:length(x)
            for jj = 1:length(y)
                if ((x(ii) <= -Cases.x_max(aa)/2) || (y(jj) == 0))
                    grad_mag(jj,ii)=0;
                end
            end
        end

    %     for ii = 1:length(x)
    %         for jj = 1:length(y)
    %             if grad_mag(jj,ii)<0.39
    %                 grad_mag(jj,ii)=0;
    %             end
    %         end
    %     end


        % Determine interface (IF) y = f(x)
        IF = zeros(3,1);
        index=2;
        % Loop all x and y
        for ii = 2:length(x)
            for jj = 2:length(y)
                if grad_mag(jj,ii)>0.9*max(max(grad_mag))
                    IF(1,index) = x(ii)+Cases.x_max(aa)/2; % Create IF profile: x
                    IF(2,index) = y(jj); % Create IF profile: y
                    index=index+1;
                end
            end
        end

        % Final deposit height
        IF(2,1)= max(IF(2,:));
        Cases.y_f_CFD(aa) = IF(2,1);

        % Final Run-out
        x_f = IF(2,:)>=Cases.d_s(aa); % Get logical array which contains zeros for all positions where the deposit height is less than the particle diameter
        % Check if logical array contains zeros
        if nnz(x_f)==length(x_f)
            Cases.x_f_CFD(aa) = IF(1,end); 
        else
            position = find(x_f==0, 1, 'first');
            if x_f(position+1) ~= 0
                position = find(x_f==0, 1, 'last');
            end
            Cases.x_f_CFD(aa) = IF(1,position);
        end

        % Plot IF on dimensional shape figure


        % Plot on dimensional shape figure
        switch aa
            case 1
                figure(fig_small);
                markersize = markersizelist(1);
                markertype = 'o';
            case 2
                figure(fig_medium);
                markersize = markersizelist(2);
                markertype = 'o';
            case 3
                figure(fig_large);
                markersize = markersizelist(3);
                markertype = 'o';
            case 4
                figure(fig_small);
                markersize = markersizelist(1);
                markertype = 'x';
            case 5
                figure(fig_medium);
                markersize = markersizelist(2);
                markertype = 'x';
            case 6
                figure(fig_large);
                markersize = markersizelist(3);
                markertype = 'x';
            otherwise  
        end
        plot(IF(1,:),IF(2,:),'-k');
        scatter(Cases.x_f_CFD(aa),0,markertype,'b');
        scatter(min(xlim),Cases.y_f_CFD(aa),markertype,'g');

        % Plot on non-dimensional f(AR) figure
        figure(fig_AR);
        x_n = (Cases.x_f_CFD(aa)-Cases.x_0(aa))/Cases.x_0(aa);
        x_n = Cases.x_f_CFD(aa)/Cases.x_0(aa);
        yyaxis left;
        scatter(Cases.AR(aa),x_n,markersize,markertype,'b');
        y_n = Cases.y_f_CFD(aa)/Cases.x_0(aa); % Normalized with x0 acc. to Lube et al. (2005)
        y_n = y_n/Cases.AR(aa); % Normalized with y0
        yyaxis right;
        scatter(Cases.AR(aa),y_n,markersize,markertype,'g');

        % Plot on non-dimensional shape figure
        % Scaling acc. to Lube et al. (2005)
        figure(fig_Shape_y_0);
        x_n = IF(1,:)/Cases.x_f_CFD(aa);
        y_n = IF(2,:)/Cases.y_f_CFD(aa);
        if aa<4
            scatter(x_n,smooth(y_n),markersize,'k'); % Plot IF
        else
            % scatter(x_n,smooth(y_n),markersize,'b'); % Plot IF
        end
    end
    
%     % Read cfd post data for final cliff-collapse
%     [x,y,vof] = importCFDpostdata([pwd '\' Cases.Name{aa} '_.csv']);
%     vof=vof';
% 
%     % Determine gradient
%     [dvofdx,dvofdy] = gradient(vof);
%     grad_mag = (dvofdx.^2+dvofdy.^2).^0.5;
% %     figure
% %     contour(x,y,grad_mag);
% %     contour(x,y,vof);
% %     hold on
% %     quiver(x,y,dvofdx,dvofdy)
% %     hold off
% 
%     % Set gradient to zero at domain boundaries
% %     for ii = 1:length(x)
% %         for jj = 1:length(y)
% %             if ((x(ii) <= -Cases.x_max(aa)/2+Cases.dx(aa)) || (y(jj) == 0))
% %                 grad_mag(jj,ii)=0;
% %             end
% %         end
% %     end
%     
%     for ii = 1:length(x)
%         for jj = 1:length(y)
%             if ((x(ii) <= -Cases.x_max(aa)/2) || (y(jj) == 0))
%                 grad_mag(jj,ii)=0;
%             end
%         end
%     end
%     
% %     for ii = 1:length(x)
% %         for jj = 1:length(y)
% %             if grad_mag(jj,ii)<0.39
% %                 grad_mag(jj,ii)=0;
% %             end
% %         end
% %     end
%     
%     
%     % Determine interface (IF) y = f(x)
%     IF = zeros(3,1);
%     index=2;
%     % Loop all x and y
%     for ii = 2:length(x)
%         for jj = 2:length(y)
%             if grad_mag(jj,ii)>0.9*max(max(grad_mag))
%                 IF(1,index) = x(ii)+Cases.x_max(aa)/2; % Create IF profile: x
%                 IF(2,index) = y(jj); % Create IF profile: y
%                 index=index+1;
%             end
%         end
%     end
%     
%     % Final deposit height
%     IF(2,1)= max(IF(2,:));
%     Cases.y_f_CFD(aa) = IF(2,1);
%     
%     % Final Run-out
%     x_f = IF(2,:)>=Cases.d_s(aa); % Get logical array which contains zeros for all positions where the deposit height is less than the particle diameter
%     % Check if logical array contains zeros
%     if nnz(x_f)==length(x_f)
%         Cases.x_f_CFD(aa) = IF(1,end); 
%     else
%         position = find(x_f==0, 1, 'first');
%         if x_f(position+1) ~= 0
%             position = find(x_f==0, 1, 'last');
%         end
%         Cases.x_f_CFD(aa) = IF(1,position);
%     end
%     
%     % Plot IF on dimensional shape figure
% 
%     
%     % Plot on dimensional shape figure
%     switch aa
%         case 1
%             figure(fig_small);
%             markersize = markersizelist(1);
%             markertype = 'o';
%         case 2
%             figure(fig_medium);
%             markersize = markersizelist(2);
%             markertype = 'o';
%         case 3
%             figure(fig_large);
%             markersize = markersizelist(3);
%             markertype = 'o';
%         case 4
%             figure(fig_small);
%             markersize = markersizelist(1);
%             markertype = 'x';
%         case 5
%             figure(fig_medium);
%             markersize = markersizelist(2);
%             markertype = 'x';
%         case 6
%             figure(fig_large);
%             markersize = markersizelist(3);
%             markertype = 'x';
%         otherwise  
%     end
%     plot(IF(1,:),IF(2,:),'-k');
%     scatter(Cases.x_f_CFD(aa),0,markertype,'b');
%     scatter(min(xlim),Cases.y_f_CFD(aa),markertype,'g');
%     
%     % Plot on non-dimensional f(AR) figure
%     figure(fig_AR);
%     x_n = (Cases.x_f_CFD(aa)-Cases.x_0(aa))/Cases.x_0(aa);
%     x_n = Cases.x_f_CFD(aa)/Cases.x_0(aa);
%     yyaxis left;
%     scatter(Cases.AR(aa),x_n,markersize,markertype,'b');
%     y_n = Cases.y_f_CFD(aa)/Cases.x_0(aa); % Normalized with x0 acc. to Lube et al. (2005)
%     y_n = y_n/Cases.AR(aa); % Normalized with y0
%     yyaxis right;
%     scatter(Cases.AR(aa),y_n,markersize,markertype,'g');
% 
%     % Plot on non-dimensional shape figure
%     % Scaling acc. to Lube et al. (2005)
%     figure(fig_Shape_y_0);
%     x_n = IF(1,:)/Cases.x_f_CFD(aa);
%     y_n = IF(2,:)/Cases.y_f_CFD(aa);
%     if aa<4
%         scatter(x_n,smooth(y_n),markersize,'k'); % Plot IF
%     else
%         % scatter(x_n,smooth(y_n),markersize,'b'); % Plot IF
%     end
%     
%     
%     
% %     % Plot on non-dimensional shape figure
% %     % Scaling acc. to de Vet et al. (2010)
% %     figure(fig_Shape_y_0);
% %     scatter([-1 x_n*Cases.x_0(aa)/Cases.y_0(aa)],[y_n*Cases.x_0(aa)/Cases.y_0(aa) 0],'r');  
% %     scatter([-1 Cases.x_f_CFD(aa)/Cases.y_0(aa)], [Cases.y_f_CFD(aa)/Cases.y_0(aa) 0],'r');
% %     plot((IF(1,:)-Cases.x_0(aa))/Cases.y_0(aa),IF(2,:)/Cases.y_0(aa),'r--'); % Plot IF
    
end



