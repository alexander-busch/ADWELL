function [  ] = CompParticleScale( d_p, H, v_set )
%CompParticleScale Summary of this function goes here
%   Detailed explanation goes here


% Parameters
global rho_f rho_s k_Cr n_Cr ny_0 ny_inf;
global lambda_Maxwell_1 lambda_Maxwell_2 lambda_3ITT_1 lambda_3ITT_2;
global fig_PS_3ITT fig_PS_Maxwell colorlist;

% Initialize v_set_old
v_set_old = zeros(length(v_set),1);

% Newtonian shear rate estimate
S = v_set ./ (0.5.*d_p);

while (max(abs(v_set-v_set_old)) >= min(v_set)/1e-3)
    % Update settling velocity
    v_set_old = v_set;

    % Newtonian shear rate estimate
    S = v_set./(0.5.*d_p);

    % Cross shear rate estimate
    c1 = -((1-n_Cr).*(ny_0.^2).*S.^(1-n_Cr))./(k_Cr.*(1+k_Cr.*S.^(1-n_Cr)).^2);
    c2 = ny_0./(1+k_Cr.*S.^(1-n_Cr));
    c3 = ny_0.*S./(1+k_Cr.*S.^(1-n_Cr));
    n_dash = S.*(c1+c2)./c3;
    S = (3.*n_dash+1)./(4.*n_dash).*S;

    % Apparent viscosity estimate based on shear rate estimate
    eta = ny_inf+(ny_0-ny_inf)./((1+k_Cr.*S).^n_Cr);

    % Particle Reynolds number estimate
    Re_p = rho_f.*v_set.*d_p./eta;
    
    % Settling velocity
    for i = 1:length(Re_p)        
        if (Re_p(i)< 0.1) % Stokes flow
            % Settling velocity (m/s)
            v_set(i) = (rho_s-rho_f)*9.81*(d_p(i))^2/18/eta(i);          
        else
            % Coefficient of Drag
            c_D = (24/Re_p(i))*(1+0.15*(Re_p(i)^0.687));

            % Settling velocity
            v_set(i) = sqrt(4*d_p(i)*9.81/3/c_D*(rho_s/rho_f-1));
        end
    end
end

% Particle Reynolds number
Re_p = rho_f.*v_set.*d_p./eta;

% Compute Power law n & K for current eta & S based on analytical
% solution
k_PL = exp(log(S).*(ny_0-ny_inf).*n_Cr.*k_Cr.*S./((ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf).*(k_Cr.*S+1))).*(ny_inf+(k_Cr.*S+1).^(-n_Cr)*ny_0-(k_Cr.*S+1).^(-n_Cr).*ny_inf);
n_PL=(-n_Cr.*k_Cr.*S.*ny_0+n_Cr.*k_Cr.*S.*ny_inf+k_Cr.*S.*(k_Cr.*S+1).^n_Cr.*ny_inf+k_Cr.*S.*ny_0-k_Cr.*S.*ny_inf+ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf)./(k_Cr.*S.*(k_Cr.*S+1).^n_Cr.*ny_inf+k_Cr.*S.*ny_0-k_Cr.*S.*ny_inf+ny_inf.*(k_Cr.*S+1).^n_Cr+ny_0-ny_inf);

% Generalized Reynolds number
Nom = rho_f.*v_set.^(2-n_PL).*d_p.^n_PL;
DeN = (24./2).^(n_PL-1) .* k_PL .* ((1+n_PL.*2)./(n_PL+n_PL.*2)).^n_PL;
Re_G = Nom ./ DeN;

% Time scales   
T_i = [kron(H',1./v_set),...
    d_p./v_set,...
    1 ./ S,...
    rho_s.*(d_p).^2./(18.*eta)];

% Deborah numbers
De_3ITT_1 = lambda_3ITT_1./T_i;
De_3ITT_2 = lambda_3ITT_2./T_i;
De_Maxwell_1 = lambda_Maxwell_1./T_i;
De_Maxwell_2 = lambda_Maxwell_2./T_i;


n=size(T_i);

for j = 1:n(2)
    
    % Fill areas, X
    X=[Re_p',fliplr(Re_p')]; % create continuous x value array for plotting
    
    figure(fig_PS_3ITT);
        Y=[De_3ITT_1(:,j)',fliplr(De_3ITT_2(:,j)')]; % create y values for out and then back
        h=fill(X,Y,colorlist{j});
        set(h,'facealpha',.5);
        set(h,'EdgeColor','none');

        % Plot Re_G
        vline(Re_p,'k:');
        
    figure(fig_PS_Maxwell);
        Y=[De_Maxwell_1(:,j)',fliplr(De_Maxwell_2(:,j)')]; % create y values for out and then back
        h=fill(X,Y,colorlist{j});
        set(h,'facealpha',.5);
        set(h,'EdgeColor','none');
    
        % Plot Re_G
        vline(Re_p,'k:');
end % of for loop

for k = 1:length(d_p)
    % Text
    txt_y = 1e-3;
    txt_txt = [num2str(d_p(k)*1000,2) ' mm'];

    figure(fig_PS_3ITT);
        text(Re_p(k),txt_y, txt_txt,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontSize',10,...
        'FontName','Arial',...
        'FontWeight','normal',...
        'Rotation',0);

    figure(fig_PS_Maxwell);
        text(Re_p(k),txt_y, txt_txt,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontSize',10,...
        'FontName','Arial',...
        'FontWeight','normal',...
        'Rotation',0);
end % of for
    
end % of function