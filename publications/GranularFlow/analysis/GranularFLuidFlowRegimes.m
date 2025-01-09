clear all;
close all;
clc;
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-12 - On the validity of the two-fluid-ktgf approach for dense gravity-driven granular flows\Figures\';

% Size for symbols representing particle diameters    
markersizelist = [5,10,15];
markersymbolslist = {'o','d','s'};
markercolorlist = {'[0 0 1]' 'none'; ... % Air
    '[0 0 1]' '[0 0 1]'; ... H2O
    '[0 1 0]' '[0 1 0]'; ... PAC3
    '[1 0 0]' '[1 0 0]'}; % PAC4

d_s = [0.0001 0.001 0.01];
rho_s = 2650;
rho_f = [1.225 998 1000 1000];

mu_0 = [0 0 7.21e-2 2.14e-1];
mu_inf = [1.79e-5 1.002e-3 1.002e-3 1.002e-3];
lambda_Cr = [0 0 1.09e-2 2.61e-2];
n_Cr = [0 0 0.586 0.608];

fig_handle = CreateFigure( 'Granular-fluid flow regimes (Bougouin and Lacaze 2018)', 'St [-]', 'r [-]', {'log' 'log'}, 'DINA5');
set(gca,'XLim',[1e-3 1e5],'YLim',[1e0 1e2]);
plot([10 10],[4 max(ylim)],'-k');
plot([10 max(xlim)],[4 4],'-k');
plot([2.5e0 10],[min(ylim) 4],'-k');
text(1e-1,1e1,'Viscous regime','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14,'FontWeight','bold','Color','k'); %'interpreter','latex'
text(1e3,2e1,'Free fall regime','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14,'FontWeight','bold','Color','k'); %'interpreter','latex'
text(1e3,2.6e0,'Inertial regime','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14,'FontWeight','bold','Color','k'); %'interpreter','latex'
text(max(xlim),4,'r_c = 4  ','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',14,'FontWeight','normal','Color','k'); %'interpreter','latex'
text(10,max(ylim),'St_c = 10  ','HorizontalAlignment','right','VerticalAlignment','bottom','Rotation',90,'FontSize',14,'FontWeight','normal','Color','k'); %'interpreter','latex'
text(10,4,'Re_c = 2.5','HorizontalAlignment','right','VerticalAlignment','bottom','Rotation',70,'FontSize',14,'FontWeight','normal','Color','k'); %'interpreter','latex'

TightFigure(gca);

for aa = 1:length(mu_0)
    % aa = 3
    % aa = 4
    for bb = 1:length(d_s)
        % bb = 2
        % bb = 3
        
        % Determine viscosity
        if aa<3
            % Newtonian, use high-shear Newtonian viscosity
            eta_f = mu_inf(aa);
        else
            %PAC, determine viscosity iteratively based on settling
            %velocity
            
            % Estimate Stokes settling velocity based on zero-shear-viscosity  [m/s]
            v_set = (rho_s-rho_f(aa))*9.81*(d_s(bb))^2/18/mu_0(aa);
            v_set_old = 0;
        
            % Accuracy
            accuracy = 1e-8;          
            
            
            % Iteratively compute settling velocity
            while abs(v_set - v_set_old) > accuracy % eps(v_set)

                % Update settling velocity
                v_set_old = v_set;
                
                % Shear rate
                SR = v_set/d_s(bb);

                % Apparent viscosity
                eta_f = mu_inf(aa)+(mu_0(aa)-mu_inf(aa))/(1+(lambda_Cr(aa)*SR)^(1-n_Cr(aa)));

                % Particle Reynolds number 
                Re_p = rho_f(aa)*v_set*d_s(bb)/eta_f;

                % Schiller & Naumann (1935) drag coefficient
                c_D = (24/Re_p)*(1+0.15*(Re_p^0.687));

                % Schiller & Naumann (1935) high Re correction
                c_D(Re_p > 1000)=0.44;

                % Settling velocity
                v_set = (4*d_s(bb)*9.81/3/c_D*(rho_s/rho_f(aa)-1))^0.5;
        
            end
            
        end
        
        
        % Stokes number
        St = (rho_s*(rho_s-rho_f(aa))*9.81*d_s(bb)^3).^0.5/(18*sqrt(2)*eta_f);
        
        % Froude number
        Fr = rho_f(aa)*9.81*d_s(bb)^3*(rho_s-rho_f(aa))/18^2/eta_f;
        %Fr2 = rho_f(aa)*(-150/3.5*0.6/(1-0.6)*eta_f/rho_f(aa) + sqrt((150/3.5*0.6/(1-0.6)*eta_f/rho_f(aa))^2 + d_s(bb)*(rho_s-rho_f(aa))*9.81))^2/9.81/d_s(bb)/(rho_s-rho_f(aa));
        
        % Fluid grain density ratio
        r = (rho_s/rho_f(aa))^0.5;
        
        plot(St,r,'o','MarkerEdgeColor',markercolorlist{aa,1},'MarkerFaceColor',markercolorlist{aa,2},'MarkerSize',markersizelist(bb));
        
        %text(St,r,num2str(Fr,2),'fontsize',6,'color','k','HorizontalAlignment','center','VerticalAlignment','top');
        %text(St,r,num2str(Fr2,'%4.4f'),'fontsize',8,'color','k','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end





% Specify Figure Size and Page Size
set(fig_handle,'PaperPositionMode','auto') %set paper pos for printing
fig_handle_pos = fig_handle.PaperPosition;
fig_handle.PaperSize = [fig_handle_pos(3) fig_handle_pos(4)];
fig_handle.PaperUnits = 'centimeters';

% Save Figure to File Format
fig_name = 'granularflowregimes';
print(fig_handle,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig_handle,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file

