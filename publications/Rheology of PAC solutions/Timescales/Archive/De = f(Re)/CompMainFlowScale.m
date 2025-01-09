function [  ] = CompMainFlowScale( d_o, d_i, d_h, L, VolFlowRate, rpm )
%CompMainFlowScale Summary of this function goes here
%   Detailed explanation goes here

% d_o=d_o(i)
% d_i=d_i(i)
% d_h=d_h(i)


% Parameters
global rho_f k_Cr n_Cr ny_0 ny_inf;
global lambda_Maxwell_1 lambda_Maxwell_2 lambda_3ITT_1 lambda_3ITT_2;
global fig_MFS_3ITT fig_MFS_Maxwell colorlist;

% x-sectional area
A = pi*((d_o*25.4/1000).^2-(d_i*25.4/1000).^2)/4;

% Bulk velocity, Annulus
U = VolFlowRate./A;

% Newtonian shear rate estimate, Annulus
S = 12 * U ./ d_h;

% Cross shear rate estimate
% c1 = -((1-n_Cr).*(ny_0.^2).*S.^(1-n_Cr))./(k_Cr.*(1+k_Cr.*S.^(1-n_Cr)).^2);
% c2 = ny_0./(1+k_Cr.*S.^(1-n_Cr));
% c3 = ny_0.*S./(1+k_Cr.*S.^(1-n_Cr));
% n_dash = S.*(c1+c2)./c3;
% S = (3.*n_dash+1)./(4.*n_dash).*S

% Apparent viscosity estimate based on shear rate estimate
eta = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*S).^n_Cr);

% Reynolds number, Annulus
Re = rho_f.*U.*d_h./eta;

% Compute Power-Law n & K for current eta & S based on Metzner &
% Reed
% tau = eta_a .* S;
% epsilon = 1e-6.*S;
% S_up = S+epsilon;
% S_down = S-epsilon;
% eta_up = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*S_up).^n_Cr);
% eta_down = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*S_down).^n_Cr);
% n_dash = (log(eta_up.*S_up)-log(eta_down.*S_down))./(log(S_up)-log(S_down))
% K_dash = tau./S.^n_dash;
% n_PL = n_dash;
% k_PL = K_dash./((3.*n_dash+1)./(4.*n_dash)).^n_dash;

% Compute Power law n & K for current eta & S based on analytical
% solution
k_PL = exp(log(S).*(ny_0-ny_inf).*n_Cr.*k_Cr.*S./((ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf).*(k_Cr.*S+1))).*(ny_inf+(k_Cr.*S+1).^(-n_Cr)*ny_0-(k_Cr.*S+1).^(-n_Cr).*ny_inf);
n_PL=(-n_Cr.*k_Cr.*S.*ny_0+n_Cr.*k_Cr.*S.*ny_inf+k_Cr.*S.*(k_Cr.*S+1).^n_Cr.*ny_inf+k_Cr.*S.*ny_0-k_Cr.*S.*ny_inf+ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf)./(k_Cr.*S.*(k_Cr.*S+1).^n_Cr.*ny_inf+k_Cr.*S.*ny_0-k_Cr.*S.*ny_inf+ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf);

% Generalized Reynolds number, Annulus
Nom = rho_f.*U.^(2-n_PL).*d_h.^n_PL;
DeN = (24./2).^(n_PL-1) .* k_PL .* ((1+n_PL.*2)./(n_PL+n_PL.*2)).^n_PL;
Re_G = Nom ./ DeN;



% Time scales   
T_i = [L ./ U,...
    d_h ./ U,...
    1 ./ (rpm/60*ones(length(Re),1)),...
    1 ./ S,...
    eta.*d_h./rho_f./U.^3];


% Deborah numbers
De_3ITT_1 = lambda_3ITT_1./T_i;
De_3ITT_2 = lambda_3ITT_2./T_i;
De_Maxwell_1 = lambda_Maxwell_1./T_i;
De_Maxwell_2 = lambda_Maxwell_2./T_i;



for j = 1:length(T_i)
%     plot(Re_G,De_3ITT_1(:,j),'color',colorlist{j});
%     plot(Re_G,De_3ITT_2(:,j),'color',colorlist{j});
    
    % Fill areas
    X=[Re',fliplr(Re')]; % create continuous x value array for plotting
    
    figure(fig_MFS_3ITT);
    Y=[De_3ITT_1(:,j)',fliplr(De_3ITT_2(:,j)')]; % create y values for out and then back
    h=fill(X,Y,colorlist{j});
    set(h,'facealpha',.5);
    set(h,'EdgeColor','none');
    
    % Plot Re_G
    vline(max(Re_G),'k:');
    
    figure(fig_MFS_Maxwell);
    Y=[De_Maxwell_1(:,j)',fliplr(De_Maxwell_2(:,j)')]; % create y values for out and then back
    h=fill(X,Y,colorlist{j});
    set(h,'facealpha',.5);
    set(h,'EdgeColor','none');
    
    % Plot Re_G
    vline(max(Re_G),'k:');
    
%     for k = 1:length(d_o)      
%         y1 = 1e-3;
%         txt1 = [num2str(d_o(k)/(25.4/1000),3)];
%         txt2 = [num2str(d_i(k)/(25.4/1000),2)];
%         txt3 = [num2str(d_h(k),1)];
%         
%         figure(fig_MFS_3ITT );
%         vline(Re_G(k),'k:');
%         text(Re_G(k),0.6*y1, txt1,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%         text(Re_G(k),y1, txt2,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%         text(Re_G(k),1.5*y1, txt3,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%              
%         figure(fig_MFS_Maxwell);
%         vline(Re_G(k),'k:');
%         text(Re_G(k),0.6*y1, txt1,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%         text(Re_G(k),y1, txt2,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%         text(Re_G(k),1.5*y1, txt3,...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle',...
%             'FontSize',8,...
%             'FontName','Arial',...
%             'FontWeight','normal',...
%             'Rotation',0);
%     end % of for loop
end % of for loop

end % of function






