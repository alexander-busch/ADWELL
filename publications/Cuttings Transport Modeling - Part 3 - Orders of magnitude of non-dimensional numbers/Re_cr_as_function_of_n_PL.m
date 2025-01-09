n = linspace(0,1);
CreateFigure( 'Critical Reynolds number as function of PL index n',...
    'Reynolds number Re_c_r [-]', 'PL index n [-]');
% Dodge & Metzner (1959)
plot(n,3250-1150*n);
% Ryan and Johnson (1959)
plot(n,6464*n./(3*n+1).^2.*(2+n).^((2+n)./(1+n)));
% Mishra and Tripathi (1975)
plot(n,2100*(4*n+2).*(5*n+3)./(3*(3*n+1).^2));
legend('Dodge & Metzner (1959)',...
    'Ryan and Johnson (1959)',...
    'Mishra and Tripathi (1975)');