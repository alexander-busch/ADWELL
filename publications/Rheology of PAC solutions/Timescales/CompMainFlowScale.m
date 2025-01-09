function [  ] = CompMainFlowScale( d_o, d_i, d_h, L, VolFlowRate, fluid, ii, color )
%CompMainFlowScale Summary of this function goes here
%   Detailed explanation goes here


% Parameters
global rho_f lambda_Cross n_Cross mu_0 mu_inf;
global lambda_el_upper lambda_el_lower lambda_th_upper lambda_th_lower;
global fig_MFS_3ITT fig_MFS_Maxwell fig_MFS_Pipkin_th fig_MFS_Pipkin_el colorlist;

% x-sectional area of annulus
A = pi*((d_o*25.4/1000).^2-(d_i*25.4/1000).^2)/4;

% Bulk velocity in annulus
U = VolFlowRate./A;

% Newtonian shear rate at the wall in annulus
SR = 12 * U ./ d_h;

% Cross shear rate estimate
% c1 = -((1-n_Cr).*(mu_0.^2).*S.^(1-n_Cr))./(k_Cr.*(1+k_Cr.*S.^(1-n_Cr)).^2);
% c2 = mu_0./(1+k_Cr.*S.^(1-n_Cr));
% c3 = mu_0.*S./(1+k_Cr.*S.^(1-n_Cr));
% n_dash = S.*(c1+c2)./c3;
% S = (3.*n_dash+1)./(4.*n_dash).*S

% Apparent viscosity estimate based on shear rate estimate
eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cross.*SR).^n_Cross);

% Reynolds number, Annulus
Re = rho_f.*U.*d_h./eta;

[ n_PL, k_PL ] = Cross2PL( lambda_Cross, n_Cross, mu_0, mu_inf, SR );

% Generalized Reynolds number, Annulus
Nom = rho_f.*U.^(2-n_PL).*d_h.^n_PL;
DeN = (24./2).^(n_PL-1) .* k_PL .* ((1+n_PL.*2)./(n_PL+n_PL.*2)).^n_PL;
Re_G = Nom ./ DeN;



% Time scales
T_MFT = L ./ U;
%T_LET = d_h ./ U;
T_SR  = 1 ./ SR;
%T_Rot = 1 ./ (rpm/60); % (rpm/60*ones(length(Re),1));
%T_Kol = eta.*d_h./rho_f./U.^3;


% Deborah numbers
figure(fig_MFS_3ITT);

X=[Re_G,fliplr(Re_G)]; % create continuous x value array for plotting

De_th_lower = lambda_th_lower./T_MFT;
De_th_upper = lambda_th_upper./T_MFT;
% plot(Re_G,De_th_lower,'color',colorlist{1});
% plot(Re_G,De_th_upper,'color',colorlist{1});
Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{1});
set(h,'facealpha',0.3);
set(h,'EdgeColor','none');


% De_th_lower = lambda_th_lower./T_LET;
% De_th_upper = lambda_th_upper./T_LET;
% % plot(Re_G,De_th_lower,'color',colorlist{2});
% % plot(Re_G,De_th_upper,'color',colorlist{2});
% Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{2});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');

De_th_lower = lambda_th_lower./T_SR;
De_th_upper = lambda_th_upper./T_SR;
% plot(Re_G,De_th_lower,'color',colorlist{3});
% plot(Re_G,De_th_upper,'color',colorlist{3});
Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{2});
set(h,'facealpha',0.3);
set(h,'EdgeColor','none');

% De_th_lower = lambda_th_lower./(max(T_Rot)*ones(length(Re),1));
% De_th_upper = lambda_th_upper./(min(T_Rot)*ones(length(Re),1));
% % plot(Re_G,De_th_lower,'color',colorlist{4});
% % plot(Re_G,De_th_upper,'color',colorlist{4});
% Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{4});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');
% 
% De_th_lower = lambda_th_lower./T_Kol;
% De_th_upper = lambda_th_upper./T_Kol;
% % plot(Re_G,De_th_lower,'color',colorlist{5});
% % plot(Re_G,De_th_upper,'color',colorlist{5});
% Y=[De_th_lower(:)',fliplr(De_th_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{5});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');

vline(max(Re_G),'k:');
y1 = 1e-2;
txt1 = num2str(ii); %d_o/(25.4/1000),3);
% txt2 = num2str(d_i/(25.4/1000),2);
% txt3 = num2str(d_h,1);
text(max(Re_G),y1*ii, txt1,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',8,...
    'FontName','Arial',...
    'FontWeight','normal',...
    'Rotation',0);
% text(max(Re_G),y1, txt2,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',8,...
%     'FontName','Arial',...
%     'FontWeight','normal',...
%     'Rotation',0);
% text(max(Re_G),1.5*y1, txt3,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',8,...
%     'FontName','Arial',...
%     'FontWeight','normal',...
%     'Rotation',0);

figure(fig_MFS_Maxwell);

De_el_lower = lambda_el_lower./T_MFT;
De_el_upper = lambda_el_upper./T_MFT;
% plot(Re_G,De_el_lower,'color',colorlist{1});
% plot(Re_G,De_el_upper,'color',colorlist{1});
Y=[De_el_lower(:)',fliplr(De_el_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{1});
set(h,'facealpha',0.3);
set(h,'EdgeColor','none');

% De_el_lower = lambda_el_lower./T_LET;
% De_el_upper = lambda_el_upper./T_LET;
% % plot(Re_G,De_el_lower,'color',colorlist{2});
% % plot(Re_G,De_el_upper,'color',colorlist{2});
% Y=[De_el_lower(:)',fliplr(De_el_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{2});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');

De_el_lower = lambda_el_lower./T_SR;
De_el_upper = lambda_el_upper./T_SR;
% plot(Re_G,De_el_lower,'color',colorlist{3});
% plot(Re_G,De_el_upper,'color',colorlist{3});
Y=[De_el_lower(:)',fliplr(De_el_upper(:)')]; % create y values for out and then back
h=fill(X,Y,colorlist{2});
set(h,'facealpha',0.3);
set(h,'EdgeColor','none');

De_el_lower/Re_G
De_el_upper/Re_G

% De_el_lower = lambda_el_lower./(max(T_Rot)*ones(length(Re),1));
% De_el_upper = lambda_el_upper./(min(T_Rot)*ones(length(Re),1));
% % plot(Re_G,De_el_lower,'color',colorlist{4});
% % plot(Re_G,De_el_upper,'color',colorlist{4});
% Y=[De_el_lower(:)',fliplr(De_el_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{4});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');
% 
% De_el_lower = lambda_el_lower./T_Kol;
% De_el_upper = lambda_el_upper./T_Kol;
% % plot(Re_G,De_el_lower,'color',colorlist{5});
% % plot(Re_G,De_el_upper,'color',colorlist{5});
% Y=[De_el_lower(:)',fliplr(De_el_upper(:)')]; % create y values for out and then back
% h=fill(X,Y,colorlist{5});
% set(h,'facealpha',0.3);
% set(h,'EdgeColor','none');

vline(max(Re_G),'k:');
y1 = 1e-3;
% txt1 = num2str(d_o/(25.4/1000),3);
% txt2 = num2str(d_i/(25.4/1000),2);
% txt3 = num2str(d_h,1);
text(max(Re_G),y1*ii, txt1,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',8,...
    'FontName','Arial',...
    'FontWeight','normal',...
    'Rotation',0);
% text(max(Re_G),y1, txt2,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',8,...
%     'FontName','Arial',...
%     'FontWeight','normal',...
%     'Rotation',0);
% text(max(Re_G),1.5*y1, txt3,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','middle',...
%     'FontSize',8,...
%     'FontName','Arial',...
%     'FontWeight','normal',...
%     'Rotation',0);



figure(fig_MFS_Pipkin_th);

De_th_lower = lambda_th_lower./T_MFT;
De_th_upper = lambda_th_upper./T_MFT;
%X=[De_th_upper,fliplr(De_th_lower)]; % create continuous x value array for plotting

Wi_th_lower = lambda_th_lower./T_SR;
Wi_th_upper = lambda_th_upper./T_SR;
%Y=[Wi_th_upper(:)',fliplr(Wi_th_lower(:)')]; % create y values for out and then back
%h=fill(X,Y,colorlist{1});
%set(h,'facealpha',0.3);
%set(h,'EdgeColor','none');
plot(De_th_lower,Wi_th_lower,'k:','linewidth',2);
plot(De_th_upper,Wi_th_upper,'k:','linewidth',2);


if fluid=='PAC2'
    xpos=max(De_th_upper);
    ypos=max(Wi_th_upper);
    align_hor='left';
    if ii<4
        align_ver='top';
    else
        align_ver='bottom';
    end
else
    xpos=min(De_th_lower);
    ypos=min(Wi_th_lower);
    align_hor='right';
    align_ver='bottom';
end

text(xpos,ypos, num2str(ii),...
    'HorizontalAlignment',align_hor,...
    'VerticalAlignment',align_ver,...
    'FontSize',10,...
    'FontName','Arial',...
    'FontWeight','normal',...
    'Rotation',0);
% Angle
% atan((max(Wi_th_upper)-min(Wi_th_upper))/(max(De_th_upper)-min(De_th_upper)))*180/pi
% Does not work, presumably because of log-log plot



figure(fig_MFS_Pipkin_el);

De_th_lower = lambda_el_lower./T_MFT;
De_th_upper = lambda_el_upper./T_MFT;
%X=[De_th_upper,fliplr(De_th_lower)]; % create continuous x value array for plotting

Wi_th_lower = lambda_el_lower./T_SR;
Wi_th_upper = lambda_el_upper./T_SR;
%Y=[Wi_th_upper(:)',fliplr(Wi_th_lower(:)')]; % create y values for out and then back
%h=fill(X,Y,colorlist{1});
%set(h,'facealpha',0.3);
%set(h,'EdgeColor','none');
plot(De_th_lower,Wi_th_lower,'k:','linewidth',2);
plot(De_th_upper,Wi_th_upper,'k:','linewidth',2);


if fluid=='PAC2'
    xpos=max(De_th_upper);
    ypos=max(Wi_th_upper);
    align_hor='left';
    if ii<4
        align_ver='top';
    else
        align_ver='bottom';
    end
else
    xpos=min(De_th_lower);
    ypos=min(Wi_th_lower);
    align_hor='right';
    align_ver='bottom';
end

text(xpos,ypos, num2str(ii),...
    'HorizontalAlignment',align_hor,...
    'VerticalAlignment',align_ver,...
    'FontSize',10,...
    'FontName','Arial',...
    'FontWeight','normal',...
    'Rotation',0);
% Angle
% atan((max(Wi_th_upper)-min(Wi_th_upper))/(max(De_th_upper)-min(De_th_upper)))*180/pi
% Does not work, presumably because of log-log plot

end % of function






