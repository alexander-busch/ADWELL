function [  ] = CompParticleScale( d_p, fluid, v_set )
%CompParticleScale Summary of this function goes here
%   Detailed explanation goes here


% Parameters
global rho_f rho_s lambda_Cross n_Cross mu_0 mu_inf FNSC_n FNSC_K;
global lambda_el_upper lambda_el_lower lambda_th_upper lambda_th_lower;
global fig_PS_3ITT fig_PS_Maxwell fig_PS_Pipkin_th fig_PS_Pipkin_el  colorlist;

% Initialize v_set_old
v_set_old = zeros(length(v_set),1);

% Newtonian shear rate estimate
SR = v_set ./ (0.5.*d_p);

while (abs(v_set-v_set_old) >= v_set*1e-3)
    
    % Update settling velocity
    v_set_old = v_set;

    % Newtonian shear rate estimate
    SR = v_set./(0.5.*d_p);

    % Cross shear rate estimate
%     c1 = -((1-n_Cr).*(mu_0.^2).*S.^(1-n_Cr))./(k_Cr.*(1+k_Cr.*S.^(1-n_Cr)).^2);
%     c2 = mu_0./(1+k_Cr.*S.^(1-n_Cr));
%     c3 = mu_0.*S./(1+k_Cr.*S.^(1-n_Cr));
%     n_dash = S.*(c1+c2)./c3;
%     S = (3.*n_dash+1)./(4.*n_dash).*S;

    % Apparent viscosity estimate based on shear rate estimate
    eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cross.*SR).^n_Cross);

    % Particle Reynolds number estimate
    Re_p = rho_f.*v_set.*d_p./eta;
    
    % Settling velocity
    for ii = 1:length(Re_p)        
        if (Re_p(ii)< 0.1) % Stokes flow
            % Settling velocity (m/s)
            v_set(ii) = (rho_s-rho_f)*9.81*(d_p(ii))^2/18/eta(ii);          
        else
            % Coefficient of Drag
            c_D = (24/Re_p(ii))*(1+0.15*(Re_p(ii)^0.687));

            % Settling velocity
            v_set(ii) = sqrt(4*d_p(ii)*9.81/3/c_D*(rho_s/rho_f-1));
        end
    end
end

% Particle Reynolds number
Re_p = rho_f.*v_set.*d_p./eta;

[ n_PL, k_PL ] = Cross2PL( lambda_Cross, n_Cross, mu_0, mu_inf, SR );


% Generalized Reynolds number
Nom = rho_f.*v_set.^(2-n_PL).*d_p.^n_PL;
DeN = (24./2).^(n_PL-1) .* k_PL .* ((1+n_PL.*2)./(n_PL+n_PL.*2)).^n_PL;
Re_G = Nom ./ DeN;

% Time scales   
T_PS = d_p./v_set;
T_SR = d_p./(2*v_set);
T_Rel = rho_s.*(d_p).^2./(18.*eta);

% Deborah numbers
figure(fig_PS_3ITT);

X=[Re_G',fliplr(Re_G')]; % create continuous x value array for plotting

De_th_lower = lambda_th_lower./T_PS;
De_th_upper = lambda_th_upper./T_PS;
% plot(Re_G,De_th_lower,'color',colorlist{1});
% plot(Re_G,De_th_upper,'color',colorlist{1});
% plot([Re_G Re_G],[De_th_lower De_th_upper],'color',colorlist{1},'linewidth',2);
Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{1});
set(h,'facealpha',1);
set(h,'EdgeColor','none');

De_th_lower = lambda_th_lower./T_SR;
De_th_upper = lambda_th_upper./T_SR;
% plot(Re_G,De_th_lower,'color',colorlist{3});
% plot(Re_G,De_th_upper,'color',colorlist{3});
% plot([Re_G Re_G],[De_th_lower De_th_upper],'color',colorlist{3},'linewidth',2)
Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{2});
set(h,'facealpha',1);
set(h,'EdgeColor','none');

De_th_lower = lambda_th_lower./T_Rel;
De_th_upper = lambda_th_upper./T_Rel;
% plot(Re_G,De_th_lower,'color',colorlist{2});
% plot(Re_G,De_th_upper,'color',colorlist{2});
% plot([Re_G Re_G],[De_th_lower De_th_upper],'color',colorlist{2},'linewidth',2)
Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{3});
set(h,'facealpha',1);
set(h,'EdgeColor','none');

for ii = 1:length(d_p)
    if ii==1 && strcmp(fluid,'PAC4')
        
    else
        vline(Re_G(ii),'k:');
        txt_y = 1e-3;
        txt_txt = [num2str(d_p(ii)*1000,2) ' mm'];
        text(Re_G(ii),txt_y, txt_txt,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontSize',10,...
            'FontName','Arial',...
            'FontWeight','normal',...
            'Rotation',0);
    end
end


figure(fig_PS_Maxwell);

De_el_upper = lambda_el_upper./T_PS;
De_el_lower = lambda_el_lower./T_PS;
% plot(Re_G,De_el_lower,'color',colorlist{1});
% plot(Re_G,De_el_upper,'color',colorlist{1});
% plot([Re_G Re_G],[De_el_lower De_el_upper],'color',colorlist{1},'linewidth',2);
Y=[De_el_upper(:)',fliplr(De_el_lower(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{1});
set(h,'facealpha',1);
set(h,'EdgeColor','none');

De_el_upper = lambda_el_upper./T_SR;
De_el_lower = lambda_el_lower./T_SR;
% plot(Re_G,De_el_lower,'color',colorlist{2});
% plot(Re_G,De_el_upper,'color',colorlist{2});
% plot([Re_G Re_G],[De_el_lower De_el_upper],'color',colorlist{3},'linewidth',2);
Y=[De_el_upper(:)',fliplr(De_el_lower(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{2});
set(h,'facealpha',1);
set(h,'EdgeColor','none');   

De_el_upper = lambda_el_upper./T_Rel;
De_el_lower = lambda_el_lower./T_Rel;
% plot(Re_G,De_el_lower,'color',colorlist{3});
% plot(Re_G,De_el_upper,'color',colorlist{3});
% plot([Re_G Re_G],[De_el_lower De_el_upper],'color',colorlist{2},'linewidth',2);
Y=[De_el_upper(:)',fliplr(De_el_lower(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{3});
set(h,'facealpha',1);
set(h,'EdgeColor','none');



for ii = 1:length(d_p)
    if ii==1 && strcmp(fluid,'PAC4')
        
    else
        vline(Re_G(ii),'k:');
        txt_y = 1e-3;
        txt_txt = [num2str(d_p(ii)*1000,2) ' mm'];
        text(Re_G(ii),txt_y, txt_txt,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontSize',10,...
            'FontName','Arial',...
            'FontWeight','normal',...
            'Rotation',0);
    end
end


figure(fig_PS_Pipkin_th);

De_th_lower = lambda_th_lower./T_PS;
De_th_upper = lambda_th_upper./T_PS;
%X=[De_th_upper,fliplr(De_th_lower)]; % create continuous x value array for plotting

Wi_th_lower = lambda_th_lower./T_SR;
Wi_th_upper = lambda_th_upper./T_SR;
%Y=[Wi_th_upper(:)',fliplr(Wi_th_lower(:)')]; % create y values for out and then back
%h=fill(X,Y,colorlist{1});
%set(h,'facealpha',0.3);
%set(h,'EdgeColor','none');
plot(De_th_lower,Wi_th_lower,'k:','linewidth',2);
plot(De_th_upper,Wi_th_upper,'k:','linewidth',2);



figure(fig_PS_Pipkin_el);

De_th_lower = lambda_el_lower./T_PS;
De_th_upper = lambda_el_upper./T_PS;
%X=[De_th_upper,fliplr(De_th_lower)]; % create continuous x value array for plotting

Wi_th_lower = lambda_el_lower./T_SR;
Wi_th_upper = lambda_el_upper./T_SR;
%Y=[Wi_th_upper(:)',fliplr(Wi_th_lower(:)')]; % create y values for out and then back
%h=fill(X,Y,colorlist{1});
%set(h,'facealpha',0.3);
%set(h,'EdgeColor','none');
plot(De_th_lower,Wi_th_lower,'k:','linewidth',2);
plot(De_th_upper,Wi_th_upper,'k:','linewidth',2);


% ?????????????????????
% lambda_FNSC = (FNSC_K./(2.*k_PL)).^(1./(FNSC_n-n_PL));
% De_el_FNSC = lambda_FNSC./T_i;
% c_d_ratio = 1-0.18.*(Re_p.*De_el_FNSC).^0.19;
% plot(Re_p,De_el_FNSC)
% yyaxis right;
% plot(Re_p,c_d_ratio)
% yyaxis left;
% figure;
% El = De_el_FNSC./Re_p;
% plot(Re_p,El);

   
end % of function