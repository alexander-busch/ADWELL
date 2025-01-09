fig = CreateFigure( '', '', '', {'lin','lin'},'DINA5' );
axis off;
fig_sub = tight_subplot(2,2,[.1 .06],[.06 .03],[.06 .02]);

% Center
c = [0 0];

% Eccentricity
e = 0;


% Geometry
d_o = 8.5; r_o = d_o/2;
d_i = 5.5; r_i = d_i/2;
d_h = d_o-d_i;

% Bed height
h_b = [0.33 0.66]*r_o;




% Loop bed heights
for ii = 1:length(h_b)
    
    % Compute effective variables for annulus and pipe for current bed
    % height
    [ w_a, w_p, h_f, A_a_f, A_p_f, d_h_a_eff, d_h_p_eff ] = CompEffFlowVar( d_i, d_o, d_h, h_b(ii) );
    
    % Plot Annulus
    if ii == 1
        subplot(fig_sub(1));
    else
        subplot(fig_sub(3));
    end
    hold on;
    axis equal;
    % Hole diameter
    rectangle('Position',[c-r_o d_o d_o],'Curvature',[1 1],...
        'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
    % Get axis limits
    limits = xlim;   
    % Bed height
    plot([-w_a/2 w_a/2],(r_o-h_b(ii)).*[-1 -1],'r','LineWidth',2);
	% Pipe diameter
    rectangle('Position',[c-r_i d_i d_i],'Curvature',[1 1],...
        'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
    
    
        
    % Plot Pipe
    if ii == 1
        subplot(fig_sub(2));
    else
        subplot(fig_sub(4));
    end
    hold on;
    axis equal;
    set(gca,'Xlim',limits);
    % Pipe diameter
    rectangle('Position',[c-d_h/2 d_h d_h],'Curvature',[1 1],...
        'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
    % Bed height
    plot([-w_p/2 w_p/2],(d_h/2-(d_h-h_f)).*[-1 -1],'r','LineWidth',2);

    
    
    % plot areas
    
end


%% Annulus --> Pipe

% fig = CreateFigure( '', '', '', {'lin','lin'},'DINA5' );
% axis off;
% fig_sub = tight_subplot(3,2,[.1 .06],[.06 .03],[.06 .02]);
% 
% % Center
% c = [0 0];
% 
% % Eccentricity
% e = 0;
% 
% % Geometry
% d_o = 8.5; r_o = d_o/2;
% d_i = 5.5; r_i = d_i/2;
% d_h = d_o-d_i;
% 
% % Bed height
% h_b = [0.01 0.33 0.66]*r_o;
% 
% 
% 
% % Loop bed heights
% for ii = 1:length(h_b)
%     
%     % Angle and outer width of annular bed
%     alpha = 2*acos(1-2.*h_b./d_o);
%     w = d_o.*sin(alpha/2);
%     
%     
%     % Plot Annulus
%     if ii == 1
%         subplot(fig_sub(1));
%     elseif ii == 2
%         subplot(fig_sub(3));
%     else
%         subplot(fig_sub(5));
%     end
%     hold on;
%     axis equal;
%     
%     % Hole diameter
%     rectangle('Position',[c-r_o 2*r_o 2*r_o],'Curvature',[1 1],...
%         'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
% 
%     % Get axis limits
%     limits = xlim;
%     
%     % Bed height
%     plot([-w(ii)/2 w(ii)/2],(r_o-h_b(ii)).*[-1 -1],'r','LineWidth',2);
%     
% 	% Pipe diameter
%     rectangle('Position',[c-r_i 2*r_i 2*r_i],'Curvature',[1 1],...
%         'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
%     axis off; 
%     
%     % Annular solids bed area A_b = f(h_b), Cayeux et al. (2014) 
%     if h_b(ii) <= (r_o-e-r_i)
%         A_b(ii) = acos((r_o-h_b(ii))/r_o)*r_o^2-(r_o-h_b(ii))*sqrt(r_o^2-(r_o-h_b(ii))^2);
%     elseif h_b(ii) <= (r_o-e+r_i)
%         A_b(ii) = acos((r_o-h_b(ii))/r_o)*r_o^2-(r_o-h_b(ii))*sqrt(r_o^2-(r_o-h_b(ii))^2) - (acos((r_o-h_b(ii)-e)/r_i)*r_i^2-(r_o-h_b(ii)-e)*sqrt(r_i^2+(r_o-h_b(ii)-e)^2));
%     else
%         A_b(ii) = acos((r_o-h_b(ii))/r_o)*r_o^2-(r_o-h_b(ii))*sqrt(r_o^2-(r_o-h_b(ii))^2)-pi*r_i^2;
%     end
%     
%     % Solids bed fraction in annulus (and pipe)
%     alpha_b = A_b/(pi/4*(d_o^2-d_i^2));
%     
%     % Hydraulic diameter of annulus = pipe diameter
%     d_h = d_o-d_i; r_h = d_h/2;
%     
%     % Solids bed area in pipe 
%     A_b_p(ii) = (pi/4*(d_h^2))*alpha_b(ii);
%     
%     % Fluids bed area in pipe 
%     A_f_p(ii) = (pi/4*(d_h^2))-A_b_p(ii);   
%     
%     % Iterate angle alpha for pipe
%     alpha = [1 1];
%     alpha_old = 0;
%     eps = 0.0001;
%     while abs(alpha(ii)-alpha_old)>eps
%         alpha_old = alpha(ii);
%         alpha(ii) = 8*A_b_p(ii)/d_h^2+sin(alpha(ii));
%         % alpha(ii)*180/pi
%     end
%        
%     % Solids bed height in pipe
%     h_b_p(ii) = r_h*(1-cos(alpha(ii)/2));
%     % alpha_p(ii) = 2*acos(1-2*h_b_p(ii)/d_h);
%     
%     % Fluids height in pipe
%     h_f = d_h-h_b_p;
%     
%     % Solids bed width in pipe
%     w_p(ii) = d_h*sin(alpha(ii)/2);
%     
%     % Plot pipe
%     if ii == 1
%         subplot(fig_sub(2));
%     elseif ii == 2
%         subplot(fig_sub(4));
%     else
%         subplot(fig_sub(6));
%     end
%     hold on;
%     axis equal;
%     set(gca,'Xlim',limits);
%     
%     % Pipe diameter
%     rectangle('Position',[c-r_h 2*r_h 2*r_h],'Curvature',[1 1],...
%         'FaceColor', 'w', 'Edgecolor','k','LineWidth',2);
%     
%     % Bed height
%     plot([-w_p(ii)/2 w_p(ii)/2],(r_h-h_b_p(ii)).*[-1 -1],'r','LineWidth',2);
%     
%     axis off; 
%     % plot area
%     
% end

