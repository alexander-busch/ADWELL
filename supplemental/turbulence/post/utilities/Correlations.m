% Correlations

if aa==1

    % Initialize legendentries
    legendentries = cell(4,1);

    % Newtonian, laminar & turbulent regime (Morrison 2013)
    Re_G = logspace(2,5);
    f_Morrison = (0.0076*(3170./Re_G).^0.165)./(1+(3170./Re_G).^7)+16./Re_G; 
    % plot(Re_G,0.9*f_Morrison,'LineStyle','--','Color',[0.6 0.6 0.6]);
    plot(Re_G,1.0*f_Morrison,'LineStyle','-','Color',[0.6 0.6 0.6]);
    % plot(Re_G,1.1*f_Morrison,'LineStyle','--','Color',[0.6 0.6 0.6]);
    legendentries{1,1} = 'f(Re), Newtonian (Morrison 2013)';

    % PL, turbulent regime (Irvine 1988)
    Re_G = logspace(3,5); Re_G(Re_G<=2000) = []; Re_G = [2050 Re_G];
    n_PL = min(cases.n_PL_scaled);
    plot(Re_G,((((2.^(n_PL+4))./(7.^(7.*n_PL))).*(4.*n_PL./(3.*n_PL+1)).^(3.*n_PL.^2))./Re_G).^(1./(3.*n_PL+1)),'LineStyle','-','Color','k','LineWidth',1);
    legendentries{2,1} = ['f_t_u_r(Re, n_P_L = ' num2str(n_PL) '), PL (Irvine 1988)'];

    % PL, turbulent regime (Dodge & Metzner 1959)
    f_turb = logspace(-3,-2); f_turb(f_turb<=3e-3)=[]; f_turb = [3e-3 f_turb];
    plot((10.^((1./sqrt(f_turb)+0.4./(n_PL.^1.2)).*(n_PL.^0.75)./4))./(f_turb.^((2-n_PL)./2)),f_turb,'LineStyle','--','Color','k','LineWidth',1);
    legendentries{3,1} = ['f_t_u_r(Re, n_P_L = ' num2str(n_PL) '), PL (Dodge & Metzner 1959)'];

    % Trans.-turb. boundary (Moody 1944)
    plot(Re_G,(1./(-1.8.*log10(((1600./Re_G)./3.7).^1.11+(6.9./Re_G))).^2)./4,'LineStyle',':','Color','k','LineWidth',1);
    legendentries{4,1} = 'Trans.-turb. boundary (Moody 1944)';

else
    
    % Initialize legendentries
    legendentries = cell(1,1);
    
    % Newtonian, laminar regime ()
    Re_G = 2*logspace(2,3);
    plot(Re_G,24./Re_G,...
        'LineStyle','-',...
        'Color','k',...
        'LineWidth',1);
    legendentries{1,1} = 'f(Re), Generalized Annular (Delplace & Leuliet 1995)';
    
end   
    
    
% Adjust axis
set(gca,...
    'xlim', [2e2 1e5],...
    'ylim', [4e-3 8e-2]);