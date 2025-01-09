function [ fig_AR_x, fig_AR_y, fig_Shape_y_0 ] = CreateScalingFigures( fluid, alpha_s0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Run-out distance and deposit height as function of Aspect Ratio
fig_AR_x  = figure('color','w','Name','x_f vs a','Units','centimeters','Position',[1 1 12 7.8]);
% DINA4 [1 1 21.0 14.8]
% DINA5 [1 1 12 8.45]
hold on;
grid on;
box on;
set(gca,...
    'XScale','lin',...
    'YScale','lin',...
    'xlim', [0.5 3.5],...
    'ylim', [0 9]);
xlabel('Aspect ratio a = y_0/x_0 [-]','FontSize',10);
ylabel('x_n_,_f = (x_f-x_0)/x_0 [-]','FontSize',10);
TightFigure(gca);


fig_AR_y = figure('color','w','Name','y_f vs a','Units','centimeters','Position',[1 1 12 7.8]);
hold on;
grid on;
box on;
set(gca,...
    'XScale','lin',...
    'YScale','lin',...
    'xlim', [0.5 3.5],...
    'ylim', [0 1]);
xlabel('Aspect ratio a = y_0/x_0 [-]','FontSize',10);
ylabel('y_n_,_f = y_f/y_0 [-]','FontSize',10);
TightFigure(gca);


if strcmp(fluid,'air')
    
    % Lube et al. (2005)
    ARrange = linspace(0,max(xlim))';
    x_n = zeros(length(ARrange),1);
    y_n = zeros(length(ARrange),1);

    for ii = 1:length(ARrange)
        % ii= 2

        % Run-out length ------------------------------------------------------

        % Smaller aspect ratios
        lambda_x = 1.6;
        n_x = 1;
        f1=lambda_x*ARrange(ii)^n_x;

        % Larger aspect ratios
        lambda_x = 2.2;
        n_x = 2/3;    
        f2=lambda_x*ARrange(ii)^n_x;

        % Blend
        loc = (1.8+2.8)/2;
        dis = 1;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        x_n(ii) = (1 - blend)*f1 + (blend*f2);


        % Height --------------------------------------------------------------

        % Smaller aspect ratios
        lambda_y = 1; % for cliff-collape
        n_y = 0 %1;
        f1=lambda_y*ARrange(ii)^n_y;
        f1 = f1/ARrange(ii); % Normalized with y0

        % Larger aspect ratios
        lambda_y = 1; % for cliff-collape
        n_y = -2/5;
        f2=lambda_y*ARrange(ii)^n_y;
        f2 = f2/ARrange(ii); % Normalized with y0

        % Blend
        loc = 1.15;
        dis = .2;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        y_n(ii) = (1 - blend)*f1 + (blend*f2);
    end

    figure(fig_AR_x);
    plot(ARrange,x_n-0.1,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n,'-','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n+0.1,'--','Color',[0.6 0.6 0.6]);

    figure(fig_AR_y);
    plot(ARrange,y_n-0.1,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n,'-','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n+0.1,'--','Color',[0.6 0.6 0.6]);

    % Bouguin and Lacaze (2018), FFR and IR
    ARrange = linspace(0,max(xlim))';
    x_n = zeros(length(ARrange),1);
    y_n = zeros(length(ARrange),1);

    for ii = 1:length(ARrange)
        % ii= 2

        % Run-out length ------------------------------------------------------

        % Smaller aspect ratios
        lambda_x = 2.7;
        dlambda = 0.3;
        n_x = 1;
        dn = 0;
        f1_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f1=lambda_x*ARrange(ii)^n_x;
        f1_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Larger aspect ratios
        lambda_x = 3.7;
        dlambda = 0.3;
        n_x = 0.64;    
        dn = 0.02;
        f2_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f2=lambda_x*ARrange(ii)^n_x;
        f2_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Blending
        loc = 2;
        dis = 0.3;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        x_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        x_n(ii) = (1 - blend)*f1 + (blend*f2);
        x_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);


        % Height --------------------------------------------------------------

        % Smaller aspect ratios
        lambda_y = 1;
        n_y = 0;
        f1=lambda_y*ARrange(ii)^n_y;
        f1 = f1/ARrange(ii); % Normalized with y0
        f1_lower = f1;
        f1_upper = f1;

        % Larger aspect ratios
        lambda_y = 0.8; 
        dlambda = 0.07;
        n_y = -0.65;
        dn = 0.04;
        f2_lower=(lambda_y-dlambda)*ARrange(ii)^(n_y-dn);
        f2_lower = f2_lower/ARrange(ii); % Normalized with y0
        f2=lambda_y*ARrange(ii)^n_y;
        f2 = f2/ARrange(ii); % Normalized with y0
        f2_upper=(lambda_y+dlambda)*ARrange(ii)^(n_y+dn);
        f2_upper = f2_upper/ARrange(ii); % Normalized with y0

        % Blending
        loc = 0.75;
        dis = 0.2;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        y_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        y_n(ii) = (1 - blend)*f1 + (blend*f2);
        y_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);

    end

    figure(fig_AR_x);
    plot(ARrange,x_n_lower,'k--');
    plot(ARrange,x_n,'k-');
    plot(ARrange,x_n_upper,'k--');

    figure(fig_AR_y);
    plot(ARrange,y_n_lower,'k--');
    plot(ARrange,y_n,'k-');
    plot(ARrange,y_n_upper,'k--');
    
elseif (strcmp(fluid,'h2o') || strcmp(fluid,'pac2'))
    
    % Bouguin and Lacaze (2018), FFR and IR
    ARrange = linspace(0,max(xlim))';
    x_n = zeros(length(ARrange),1);
    y_n = zeros(length(ARrange),1);

    for ii = 1:length(ARrange)
        % ii= 2

        % Run-out length ------------------------------------------------------

        % Smaller aspect ratios
        lambda_x = 2.7;
        dlambda = 0.3;
        n_x = 1;
        dn = 0;
        f1_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f1=lambda_x*ARrange(ii)^n_x;
        f1_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Larger aspect ratios
        lambda_x = 3.7;
        dlambda = 0.3;
        n_x = 0.64;    
        dn = 0.02;
        f2_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f2=lambda_x*ARrange(ii)^n_x;
        f2_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Blending
        loc = 2;
        dis = 0.3;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        x_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        x_n(ii) = (1 - blend)*f1 + (blend*f2);
        x_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);


        % Height --------------------------------------------------------------

        % Smaller aspect ratios
        lambda_y = 1;
        n_y = 0;
        f1=lambda_y*ARrange(ii)^n_y;
        f1 = f1/ARrange(ii); % Normalized with y0
        f1_lower = f1;
        f1_upper = f1;

        % Larger aspect ratios
        lambda_y = 0.8; 
        dlambda = 0.07;
        n_y = -0.65;
        dn = 0.04;
        f2_lower=(lambda_y-dlambda)*ARrange(ii)^(n_y-dn);
        f2_lower = f2_lower/ARrange(ii); % Normalized with y0
        f2=lambda_y*ARrange(ii)^n_y;
        f2 = f2/ARrange(ii); % Normalized with y0
        f2_upper=(lambda_y+dlambda)*ARrange(ii)^(n_y+dn);
        f2_upper = f2_upper/ARrange(ii); % Normalized with y0

        % Blending
        loc = 0.75;
        dis = 0.2;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        y_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        y_n(ii) = (1 - blend)*f1 + (blend*f2);
        y_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);

    end

    figure(fig_AR_x);
    plot(ARrange,x_n_lower,'k--');
    plot(ARrange,x_n,'k-');
    plot(ARrange,x_n_upper,'k--');

    figure(fig_AR_y);
    plot(ARrange,y_n_lower,'k--');
    plot(ARrange,y_n,'k-');
    plot(ARrange,y_n_upper,'k--');    
    
       
    % Bouguin and Lacaze (2018), VR
    ARrange = linspace(0,max(xlim))';
    x_n = zeros(length(ARrange),1);
    y_n = zeros(length(ARrange),1);

    for ii = 1:length(ARrange)
        % ii= 2

        % Run-out length ------------------------------------------------------

        % Smaller aspect ratios
        lambda_x = 1.5;
        dlambda = 0.1;
        n_x = 1;
        dn = 0;
        f1_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f1=lambda_x*ARrange(ii)^n_x;
        f1_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Larger aspect ratios
        lambda_x = 1.9;
        dlambda = 0.1;
        n_x = 0.64;    
        dn = 0.02;
        f2_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f2=lambda_x*ARrange(ii)^n_x;
        f2_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Blending
        loc = 2;
        dis = 0.3;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        x_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        x_n(ii) = (1 - blend)*f1 + (blend*f2);
        x_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);


        % Height --------------------------------------------------------------

        % Smaller aspect ratios
        lambda_y = 1;
        n_y = 0;
        f1=lambda_y*ARrange(ii)^n_y;
        f1 = f1/ARrange(ii); % Normalized with y0
        f1_lower = f1;
        f1_upper = f1;

        % Larger aspect ratios
        lambda_y = 0.87; 
        dlambda = 0.03;
        n_y = -0.52;
        dn = 0.02;
        f2_lower=(lambda_y-dlambda)*ARrange(ii)^(n_y-dn);
        f2_lower = f2_lower/ARrange(ii); % Normalized with y0
        f2=lambda_y*ARrange(ii)^n_y;
        f2 = f2/ARrange(ii); % Normalized with y0
        f2_upper=(lambda_y+dlambda)*ARrange(ii)^(n_y+dn);
        f2_upper = f2_upper/ARrange(ii); % Normalized with y0

        % Blending
        loc = 0.75;
        dis = 0.2;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        y_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        y_n(ii) = (1 - blend)*f1 + (blend*f2);
        y_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);

    end

    figure(fig_AR_x);
    plot(ARrange,x_n_lower,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n_upper,'--','Color',[0.6 0.6 0.6]);

    figure(fig_AR_y);
    plot(ARrange,y_n_lower,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n,'-','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n_upper,'--','Color',[0.6 0.6 0.6]);

    
else
        
    % Bouguin and Lacaze (2018), VR
    ARrange = linspace(0,max(xlim))';
    x_n = zeros(length(ARrange),1);
    y_n = zeros(length(ARrange),1);

    for ii = 1:length(ARrange)
        % ii= 2

        % Run-out length ------------------------------------------------------

        % Smaller aspect ratios
        lambda_x = 1.5;
        dlambda = 0.1;
        n_x = 1;
        dn = 0;
        f1_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f1=lambda_x*ARrange(ii)^n_x;
        f1_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Larger aspect ratios
        lambda_x = 1.9;
        dlambda = 0.1;
        n_x = 0.64;    
        dn = 0.02;
        f2_lower=(lambda_x-dlambda)*ARrange(ii)^(n_x-dn);
        f2=lambda_x*ARrange(ii)^n_x;
        f2_upper=(lambda_x+dlambda)*ARrange(ii)^(n_x+dn);

        % Blending
        loc = 2;
        dis = 0.3;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        x_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        x_n(ii) = (1 - blend)*f1 + (blend*f2);
        x_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);


        % Height --------------------------------------------------------------

        % Smaller aspect ratios
        lambda_y = 1;
        n_y = 0;
        f1=lambda_y*ARrange(ii)^n_y;
        f1 = f1/ARrange(ii); % Normalized with y0
        f1_lower = f1;
        f1_upper = f1;

        % Larger aspect ratios
        lambda_y = 0.87; 
        dlambda = 0.03;
        n_y = -0.52;
        dn = 0.02;
        f2_lower=(lambda_y-dlambda)*ARrange(ii)^(n_y-dn);
        f2_lower = f2_lower/ARrange(ii); % Normalized with y0
        f2=lambda_y*ARrange(ii)^n_y;
        f2 = f2/ARrange(ii); % Normalized with y0
        f2_upper=(lambda_y+dlambda)*ARrange(ii)^(n_y+dn);
        f2_upper = f2_upper/ARrange(ii); % Normalized with y0

        % Blending
        loc = 0.75;
        dis = 0.2;
        blend = (1+ tanh((ARrange(ii)-loc)/dis))/2;
        y_n_lower(ii) = (1 - blend)*f1_lower + (blend*f2_lower);
        y_n(ii) = (1 - blend)*f1 + (blend*f2);
        y_n_upper(ii) = (1 - blend)*f1_upper + (blend*f2_upper);

    end

    figure(fig_AR_x);
    plot(ARrange,x_n_lower,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,x_n_upper,'--','Color',[0.6 0.6 0.6]);

    figure(fig_AR_y);
    plot(ARrange,y_n_lower,'--','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n,'-','Color',[0.6 0.6 0.6]);
    plot(ARrange,y_n_upper,'--','Color',[0.6 0.6 0.6]);

end



% figure(fig_AR);
% legend('x_n (Lube et al. 2005)',...
%     'y_n (Lube et al. 2005)',...
%     'Location','northwest');


%% Deposit height as function of Run-out distance, i.e. shape of final profiles

fig_Shape_y_0  = figure('color','w','Name','y_n vs x_n','Units','centimeters','Position',[10 10 12 8.45]);
hold on;
grid on;
box on;
axis equal;
set(gca,...
    'XScale','lin',...
    'YScale','lin',...
    'xlim', [0 1],...
    'ylim', [0 1],...
    'FontSize',10);

xlabel('Run-out distance x_n = x/x_0 [-]');
ylabel('Deposit height y_n = y/y_0 [-]');

figure(fig_Shape_y_0);


%% Final shapes of Lube et al. (2005)

% Import data
[~, ~, raw] = xlsread([pwd '\Shape_Lube2005.xlsx'],'Sheet1','A3:B82');
data = reshape([raw{:}],size(raw));
x_n = data(:,1);
y_n = data(:,2);
clearvars data raw;
plot(x_n,y_n);


%% Final shapes of de Vet et al. (2010)

% y_n = linspace(0,1);
% xlabel('Run-out distance x_n = (x-x_0)/y_0 [-]');
% ylabel('Deposit height y_n = y/y_0 [-]');


% Definitions
% AR = y_0/x_0; % Aspect ratio
% x_n = x/y_0; % Dimensionless x
% y_n = y/y_0; % Dimensionless y
% t_n = t/sqrt(y_0/g);

% Fitted coefficients based on fit function of Savoshi and Kudrolli (2005)
% a = -1.992;
% b = 0.050;
% x0 = 0.944;
% x_n = a.*y_n-b.*log(y_n)+x0;
% plot(x_n,y_n);

% Fitted coefficients based on fit function of de Vet and de Bruyn (2007)
% yc = -0.20;
% c = 0.511;
% x0 =1.28;
% y_n = yc+sqrt(yc.^2+c.^2*(x_n-x0).^2);
% plot(x_n,y_n);

%% Final shapes of Siavoshi and Kudrolli (2005)
% figure(fig_Shape_y_0);
% 
% % Fitted coefficients based on fit function of Savoshi and Kudrolli (2005)
% a = -1.865;
% b = 0.05;
% x0 = 1.877;
% x_n = a.*y_n-b.*log(y_n)+(x0-1);
% plot(x_n,y_n);


%% Experiment 1 of Mangeney et al. (2010)
% d_s = 600e-6; % 800e-6
% rho_s = 2500; % kg/m³
% theta_r = 23.5; % repose angle
% theta_a = 25.5; % avalanche angle
% x_0 = 0.20; % Initial width
% y_0 = 0.14; % Initial height
% AR = y_0/x_0; % Aspect ratio
% alpha_s_0 = 0.62; % Initial volume fraction
% 
% % Scaling law coefficients
% k = 0.5; % 0.4
% beta = 0.76;
% %n_x = 1.0; % For AR < 1.6 acc. to Delannay et al. (2017)
% n_y = 2/3;
% 
% % Dimensional run-out distance
% x_f = y_0*2*k/tan(theta_a*pi/180); % Equation (5)
% % Dimensional deposit height
% y_f = y_0*beta*(tan(theta_a*pi/180)/AR/k)^n_y; % Equation (6)
% 
% % Normalized run-out distance
% x_n = (x_f-x_0)/x_0; % Coordinate system with x0 = 0
% x_n = x_f/x_0; % Coordinate system with x0 = x_0
% figure(fig_AR);
% scatter(AR,x_n,'b');
% % Normalized deposit height
% y_n = y_f/x_0;
% figure(fig_AR);
% scatter(AR,y_n,'g');
% 
% figure(fig_Shape_y_0);
% scatter([-1 x_n*x_0/y_0-1],[y_n*x_0/y_0 0],'r');
% 
% % Dimensional final run-out distance x_f and final height y_f based on
% % experimental result
% x_f = 0.275;
% x_y = 0.14;
% 
% % Normalized run-out distance
% x_n = (x_f-x_0)/x_0; % Coordinate system with x0 = 0
% x_n = x_f/x_0; % Coordinate system with x0 = x_0
% figure(fig_AR);
% scatter(AR,x_n,'b+');
% % Normalized deposit height
% y_n = y_f/x_0;
% figure(fig_AR);
% scatter(AR,y_n,'g+');
% 
% figure(fig_Shape_y_0);
% scatter([-1 x_n*x_0/y_0-1],[y_n*x_0/y_0 0],'r+');

end
