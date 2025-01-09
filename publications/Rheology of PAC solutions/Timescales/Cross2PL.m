function [ n_PL, K_PL ] = Cross2PL( lambda_Cr, n_Cr, mu_0, mu_inf, SR )
%   Cross2PL Compute PL coefficients from Cross rheological model
%   !!! Only valid for Cross model of the form eta = mu_inf+(mu_0-mu_inf)./(1+(lambda_Cr.*SR).^n_Cr);

% Solution for n_PL & K_PL > 0
K_PL = (SR.^(n_Cr.*(mu_inf.*lambda_Cr.^(2.*n_Cr).*SR.^(2.*n_Cr)+2.*mu_0.*lambda_Cr.^n_Cr.*SR.^n_Cr+mu_0)./(mu_inf.*lambda_Cr.^(2.*n_Cr).*SR.^(2.*n_Cr)+SR.^n_Cr.*(mu_0+mu_inf).*lambda_Cr.^n_Cr+mu_0)).*mu_inf.*lambda_Cr.^n_Cr+SR.^(lambda_Cr.^n_Cr.*SR.^n_Cr.*n_Cr.*(mu_0-mu_inf)./(mu_inf.*lambda_Cr.^(2.*n_Cr).*SR.^(2.*n_Cr)+SR.^n_Cr.*(mu_0+mu_inf).*lambda_Cr.^n_Cr+mu_0)).*mu_0)./(1+lambda_Cr.^n_Cr.*SR.^n_Cr);
n_PL = (mu_inf.*lambda_Cr.^(2.*n_Cr).*SR.^(2.*n_Cr)-((-n_Cr-1).*mu_inf+mu_0.*(n_Cr-1)).*SR.^n_Cr.*lambda_Cr.^n_Cr+mu_0)./(mu_inf.*lambda_Cr.^(2.*n_Cr).*SR.^(2.*n_Cr)+SR.^n_Cr.*(mu_0+mu_inf).*lambda_Cr.^n_Cr+mu_0);

end



