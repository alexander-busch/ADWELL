clc;
close all;
clear all;
addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange';

filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive\180718_1';

% Aspect ratio
a = 2;
% a = 3;

% Initial solid volume fraction
alpha_s0 = 0.55;


fluid = 'pac4'; % air water pac2 pac4
num_exp = 10; % Number of exported data files per case, air 1-4, liquids 1-10



fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-12 - On the validity of the two-fluid-ktgf approach for dense gravity-driven granular flows\Figures\';

fig_title = ['Dimensional height y as function of run-out distance x and time t, sand in ' fluid ];
fig_handle = CreateFigure( fig_title, '', '', 'off', [1 1 29.7 15]);
fig_sub = tight_subplot(3,3,[.1 .06],[.06 .03],[.06 .02]);



% Size for symbols representing particle diameters    
markersizelist = [30,80,200];
markersymbolslist = {'o','d','s'};

CreateScalingFigures(fluid, alpha_s0 );

% Create case table
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
Cases.alpha_s_0 = 0.61*ones(length(index),1);   % Initial average volume fraction


% Loop all cases
for aa = 1:length(index)
    % Debugging
    % aa = 2
    % aa = 3

    % Create dimensional deposit shape figure or select correct figure
    %fig_name = ['Case ' num2str(aa) ' - Dimensional deposit height as function of run-out distance (Shape of final profiles)'];
    %fig_handle  = figure('color','w','Name',fig_name,'Units','centimeters','Position',[10 10 12 8.45]);
    figure(fig_handle);
    subplot(fig_sub(aa));
    hold on;
    plot([0 Cases.x_0(aa)],[Cases.y_0(aa) Cases.y_0(aa)],'k--');
    plot([Cases.x_0(aa) Cases.x_0(aa)],[0 Cases.y_0(aa)],'k--');
    axis equal;
    box on;
    grid on;
    set(gca,...
        'XScale','lin',...
        'YScale','lin',...
        'xlim', [0 Cases.x_max(aa)],...
        'ylim', [0 Cases.y_max(aa)]);

    fig_vof = figure; axis equal; hold on;
    set(gca,'XLim',[0 Cases.x_max(aa)],'YLim',[0 Cases.y_max(aa)]);
	fig_grad = figure; axis equal; hold on;
    set(gca,'XLim',[0 Cases.x_max(aa)],'YLim',[0 Cases.y_max(aa)]);
    
    for bb = 0:num_exp
        % Change length of loop depending 
        % Debugging
        % bb = 1
        % bb = 2
        % bb = 3
        % bb = 4
        % bb = 10

        % Variant 1 - Read Fluent exports, unfinished
        
%         % Read Fluent data
%         [x,y,vof] = importFLUENTdata([filepath '\' num2str(aa) '\solution_case\exports\vof_final_cc.txt']);
%         
%         % Shift x-data such that x0 = 0
%         x = x + Cases.x_max(aa)/2;
%                 
%         % Fit 2D surface
%         [fitresult, gof] = create3Dfit(x, y, vof, 'cubicinterp');
%         
%         % Gradient
%         [fx, fy] = differentiate(fitresult, x, y)
%         
%         % Gradient magnitude
%         grad_mag = (fx.^2+fy.^2).^0.5;
%         
%         % Fit surface to gradient magnitude
%         [fitresult, gof] = create3Dfit(x, y, grad_mag, 'cubicinterp');
        
        
        
        % Variant 2 - Read CFDpost exports, works
                
        % Read cfd post data
        if num_exp<10
            [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports\a=' num2str(a) '\' fluid '\case' num2str(Cases.Index(aa)) '_t' num2str(bb) '.csv']);     
        else
            [x,y,vof] = importCFDpostdata([pwd '\cfdpostexports\a=' num2str(a) '\' fluid '\case' num2str(Cases.Index(aa)) '_t' num2str(bb,'%02d') '.csv']);     
        end

        % Shift x-data such that x0 = 0
        x = x + Cases.x_max(aa)/2;

        % Transpose vof and remove zeros
        vof=vof';
        vof(vof==0)=nan;

        % Fill missing data points, does work but does not help 
        % vof = inpaintn(vof)
        
        %Tried but does not work
        %vof=interp2(vof,X,Y);
        %vof=fillmissing(vof);
        %dvofdy=smooth3(dvofdy)
        
        figure(fig_vof); hold on;
        % surf(x,y,vof);
        [X,Y] = meshgrid(x,y);
        plot3(X,Y,vof,'.');
                
        % Determine gradient
        [dvofdx,dvofdy] = gradient(vof);
        grad_mag = (dvofdx.^2+dvofdy.^2).^0.5;

        figure(fig_grad); hold on;
        [X,Y] = meshgrid(x,y);
        plot3(X,Y,grad_mag,'.');
        
%         figure;
%         contour(x,y,grad_mag);
%         contour(x,y,vof);
%         hold on
%         quiver(x,y,dvofdx,dvofdy)
%         hold off

        % Determine interface (IF) y = f(x)
        IF = zeros(length(x),3);
        % Loop all x and evaluate y-direction vectors 
        for ii = 1:length(x)
            %ii=65
            IF(ii,1) = x(ii); %+Cases.x_max(aa)/2; % Create IF profile: x
           
            % Current vof vector in y-direction
            current_vof = vof(:,ii);
            
            if ((max(current_vof)<0.5) || (isnan(max(current_vof))))
                IF(ii,2) = nan;
            else
                IF(ii,2) = max(y(current_vof>0.5));
            end
            
            % Current grad_mag vector in y-direction
            current_grad_mag = grad_mag(:,ii);
            if isnan(max(current_grad_mag))
                IF(ii,3) = nan;
            else
                IF(ii,3) = y(current_grad_mag==max(current_grad_mag)); % Create IF profile: y      
            end
        end
        
        % Indices of unique values in IF(:,1) = x    
        [~, ind] = unique(IF(:,1), 'rows');
        % Interpolate IF in order to close the data gap @ x = 0.05, 0.05, 5
        xi = linspace(0,Cases.x_max(aa));
        yi = interp1(IF(ind,1),IF(ind,3),xi);
        
        % Number of cells for width of initial conditions case
        cells = round(Cases.x_0(aa)/Cases.dx(aa));
        % Final deposit height
        Cases.y_f_CFD(aa) = max(yi(1:cells));
                
        % Remove wrong data points at end
        logic = yi>Cases.y_f_CFD(aa);
        yi(logic) = nan;
               
        % Plot on VOF figure
        figure(fig_vof);
        plot3(IF(:,1),IF(:,2),0.65*ones(length(x),1),'linewidth',3);
        
        % Plot on gradVOF figure
        figure(fig_grad);
        plot3(IF(:,1),IF(:,3),0.2*ones(length(x),1),'linewidth',3);
        plot3(xi,yi,0.21*ones(length(xi),1),'linewidth',3);
                        
        % Final run-out
        if aa~=3
            logic = yi>=Cases.d_s(aa); % Get logical array which contains zeros for all positions where the deposit height is less than the particle diameter
            % Check if logical array contains zeros
            if nnz(logic)==length(logic)
                Cases.x_f_CFD(aa) = xi(end); 
            else
                position = find(logic==0, 1, 'first');           
                if logic(position+1) ~= 0
                    position = find(logic==0, 1, 'last');
                end
                Cases.x_f_CFD(aa) = xi(position);
            end
        else
            position = find(isnan(yi), 1, 'first')-1;
            Cases.x_f_CFD(aa) = xi(position);
        end
        
        % Plot on other figures
        figure(fig_handle);
        subplot(fig_sub(aa));

        if bb==0

            % Effective y at t=0 
            Cases.y_0(aa) = yi(1);

            % Plot IC
            txt='k-';
            plot([0 Cases.x_0(aa)],[Cases.y_0(aa) Cases.y_0(aa)],txt); % Plot IF
            plot([Cases.x_0(aa) Cases.x_0(aa)],[0 Cases.y_0(aa)],txt); % Plot IF

            % Effective aspect ratio for settled column
            Cases.AR(aa) = Cases.y_0(aa)/Cases.x_0(aa);
                
        elseif bb==num_exp
            
            % Plot on dimensional shape figure
            plot(xi(:),yi(:),'-k');               

            % Select marker type and color
            if aa<4
                markersymbol = markersymbolslist{1};
            elseif aa<7
                markersymbol = markersymbolslist{2}; 
            else
                markersymbol = markersymbolslist{3};
            end

            % Select marker size
            if Cases.d_s(aa)==Cases.d_s(1)
                markersize = markersizelist(1);
            elseif Cases.d_s(aa)==Cases.d_s(2)
                markersize = markersizelist(2);
            else
                markersize = markersizelist(3);
            end

            % Mark x_f
            scatter(Cases.x_f_CFD(aa),0,markersize,markersymbol,'b');
            % Mark y_f
            scatter(min(xlim),Cases.y_f_CFD(aa),markersize,markersymbol,'g');

            % Plot on non-dimensional shape figure
            % Scaling acc. to Lube et al. (2005)
            figure(fig_Shape_y_0);
            x_n = xi(:)/Cases.x_f_CFD(aa);
            y_n = yi(:)/Cases.y_f_CFD(aa);
            scatter(x_n,smooth(y_n),markersymbol,'k'); % Plot IF

            % Plot on non-dimensional f(AR) figure
            figure(fig_AR);
            x_n = (Cases.x_f_CFD(aa)-Cases.x_0(aa))/Cases.x_0(aa);
            x_n = Cases.x_f_CFD(aa)/Cases.x_0(aa);
            yyaxis left;
            scatter(Cases.AR(aa),x_n,markersize,markersymbol,'b');
            y_n = Cases.y_f_CFD(aa)/Cases.x_0(aa); % Normalized with x0 acc. to Lube et al. (2005)
            y_n = y_n/Cases.AR(aa); % Normalized with y0
            yyaxis right;
            scatter(Cases.AR(aa),y_n,markersize,markersymbol,'g');

        else
            % Plot IF on dimensional shape figure           
            plot(xi(:),yi(:),'-k');
        end
    end
    
    close(fig_vof);
    close(fig_grad);   

    
% yyaxis left;
%     set(gca,...
%     'XScale','lin',...
%     'YScale','lin',...
%     'xlim', [0 4],...
%     'ylim', [0 10]);
    
    
end


% Export AR figure
figure(fig_AR);

% Specify Figure Size and Page Size
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = ['a=' num2str(a) '_' fluid '_scaling'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file


% Export dimensional shape figure 
figure(fig_handle);

% Specify Figure Size and Page Size
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = ['a=' num2str(a) '_' fluid '_shape_dimensional'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

% Export dimensional shape figure 
figure(fig_Shape_y_0);


% Specify Figure Size and Page Size
fig=gcf;
set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% fig.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = ['a=' num2str(a) '_' fluid '_shape_non-dimensional'];
print(fig,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file