
close all;
clear all;
clc;


do = 0.08;
di = 0.075;
rho = 1000; 
b = 3; % Pipe 3, Annulus 2
mdot = linspace(0.1,1);
n = [1; 0.86; 0.8; 0.76];
K = [0.00102; 0.025; 0.072; 0.118];

colorlist = {'[0 0 1]' ... % Blue 
    '[0 1 1]' ... % Cyan [0 1 1]
    '[0 1 1]' ... % Cyan [0 1 1]
    '[0 1 1]' ... % Cyan [0 1 1]
    '[1 0 1]' ... % Magenta
    '[1 0 0]'}; % Red


figure('color','w');

hold on;

for i=1:length(n)
    Usl = SuperficialVelocity(mdot, rho, do, di, 2);
    ReG_annulus = GeneralizedReynoldsNumber(rho, Usl, do, di, 2, n(i), K(i));
    
    if n==1
        L_e = 1.359.*(do-di).*ReG_annulus.^0.25; % H2O, turb, https://en.wikipedia.org/wiki/Entrance_length
    else
        L_e = (do-di).*((0.25.*n(i).^2-0.675.*n(i)+1.03).^1.6+(0.0567.*ReG_annulus).^1.6).^(1./1.6); % PAC, lam, Chhabra (2011), page 189
    end
      
    Usl = SuperficialVelocity(mdot, rho, do, di, 3);
    ReG_pipe = GeneralizedReynoldsNumber(rho, Usl, do, di, 3, n(i), K(i));
    
    subplot(1,2,1);
    hold on;
    % yyaxis left
    plot(mdot,ReG_annulus);
    % yyaxis right
    % plot(mdot,L_e);
    subplot(1,2,2);
    hold on;
    plot(mdot,ReG_pipe);
    
    legendInfo{i} = ['PAC' num2str(i)-1 ' (n = ' num2str(n(i)) ', K = '  num2str(K(i)) ')'];
end

legendInfo{1} = ['H2O (n = ' num2str(n(1)) ', K = '  num2str(K(1)) ')'];
titlename = {'Annular section', 'Pipe section'};

for i=1:2
    subplot(1,2,i);
    hold on;
    grid('on');
    % yyaxis left
    set(gca,...
        'box','on',...
        'FontSize',18,...
        'xlim', [0.2 1],...
        'ylim', [1e1 1e5],...
        'XScale','log',...
        'YScale','log',...
        'XTick',[0.2 0.5 1],...
        'YTick',[1e1 1e2 1e3 1e4 1e5]);
    xlabel('Fluid mass flow rate [kg/s]');
    ylabel('Generalized Reynolds number [-]');
    legend(legendInfo);
    title(titlename(i));

end

