i=5
d_o = d_o(i)
d_i = d_i(i)
d_h = d_h(i)

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



% Time scale
% T_i = [L ./ U,...
%     d_h ./ U,...
%     1 ./ (rpm/60*ones(length(Re),1)),...
%     1 ./ SR,...
%     eta.*d_h./rho_f./U.^3];
T_i = L./U;


% Deborah numbers
De_th_lower = lambda_th_lower./T_i;
De_th_upper = lambda_th_upper./T_i;
De_el_upper = lambda_el_upper./T_i;
De_el_lower = lambda_el_lower./T_i;


f_T_i = logspace(-1,4);


for i = 1:length(f_T_i)
    De_th_lower(i,:) = lambda_th_lower./(f_T_i(i).*T_i);
    De_th_upper(i,:) = lambda_th_upper./(f_T_i(i).*T_i);
end


fig_name = 'Annulus scale, Thixotropy';

fig_color = 'w';
fig_units = 'centimeters';
fig_position = [1 1 21 14.8]; % DIN A5 

fig_xlabel = 'Re [-]';
fig_ylabel = 'T_o_b_s/T [-]';
fig_zlabel = 'De_t_h [-]';

fig_xlim = [min(Re) max(Re)];
fig_ylim = [min(f_T_i) max(f_T_i)];
fig_zlim = [0.1 max(max(max(De_th_lower)),max(max(De_th_upper)))];

fig_scale = 'log';

Create_Figure3D(fig_name, fig_color, fig_units, fig_position, fig_xlabel, fig_ylabel, fig_zlabel, fig_xlim, fig_ylim, fig_zlim, fig_scale);
surf(Re,f_T_i,De_th_lower,'FaceColor','[1 .5 0]');
surf(Re,f_T_i,De_th_upper,'FaceColor','red');
view([-1,-1,1])
