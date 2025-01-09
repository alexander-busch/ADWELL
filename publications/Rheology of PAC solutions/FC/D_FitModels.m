clc;

%% Create figure

close all;
fig = figure; % ('units','normalized','outerposition',[0 0 1 1]);
hold on;

%% Formating

xlabel('Shear rate [1/s]');
ylabel('Apparent viscosity [Pa.s]');
grid('on');
set(gca,...
    'XScale','log',...
    'YScale','log',...
    'xlim', [0.01 1200],...
    'ylim', [9e-3 2e-0],... %3e-1
    'XTick', [1e-2 1e-1 1e-0 1e1 1e2 1e3],...
    'YTick', [1e-2 2e-2 3e-2 4e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1e-0 2e-0],...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');


%% PAC2

% PAC 2 - Laboratory 2 - UiS - insufficient shear rate range
addpath([pwd '\Data\Lab2_UiS'])
load pac2_1;
load pac2_2;
pac2(1,:) = pac2_1(:,2);
pac2(2,:) = pac2_2(:,2);
%pac2 = mean(pac2,2);
pac2 = (pac2(1,:)+pac2(2,:))./2;
%Lab2_PAC2_M_Up = plot(pac2_1(1:25),pac2(1:25)','-','Color',col_PAC2);
%Lab2_PAC2_M_Down = plot(pac2_1(26:50),pac2(26:50)','-','Color',col_PAC2);
pac2_mean=(pac2(1:25)'+flip(pac2(26:50)'))./2;
plot(pac2_1(1:25),pac2_mean,'o','Color',col_PAC2,'markersize',MS,'MarkerFaceColor','w');


%% Fit models
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\viscosity-based'

% Cross
[Cross Cross_gof] = createFit_Cross(pac2_1(1:25),pac2_mean);
close gcf;
plot(SR,Cross.mu_inf+(Cross.mu_0-Cross.mu_inf)./(1+(Cross.lambda.*SR).^Cross.n));
mu_0_PAC2 = Cross.mu_0;
mu_inf_PAC2 = Cross.mu_inf;
lambda_Cr_PAC2 = Cross.lambda;
n_Cr_PAC2 = Cross.n;

% Carreau
[Carreau Carreau_gof] = createFit_Carreau(pac2_1(1:25),pac2_mean);
close gcf;
plot(SR,Carreau.mu_inf+(Carreau.mu_0-Carreau.mu_inf).*(1+(Carreau.lambda.*SR).^2).^((Carreau.n-1)./2));
lambda_Ca_PAC2 = Carreau.lambda;


% PL - Time constant version
[PL PL_gof] = createFit_PowerLaw_TimeConstantVersion(pac2_1(1:25),pac2_mean,Carreau.mu_0);
close gcf;
plot(SR,PL.mu0.*(PL.lambda.*SR).^(PL.n-1));

[PL_lowSR PL_lowSR_gof] = createFit_PowerLaw_TimeConstantVersion(pac2_1(1:7),pac2_mean(1:7),Carreau.mu_0);
close gcf;
plot(SR,PL_lowSR.mu0.*(PL_lowSR.lambda.*SR).^(PL_lowSR.n-1));

[PL_highSR PL_highSR_gof] = createFit_PowerLaw_TimeConstantVersion(pac2_1(17:25),pac2_mean(17:25),Carreau.mu_0);
close gcf;
plot(SR,PL_highSR.mu0.*(PL_highSR.lambda.*SR).^(PL_highSR.n-1));

% Time scales
lambda_PL_PAC2 = PL.lambda;
lambda_PL_lowSR_PAC2 = PL_lowSR.lambda;
lambda_PL_highSR_PAC2 = PL_highSR.lambda;
    
% Ordinary PL coefficients
k_PAC2 = PL.mu0*PL.lambda^(PL.n-1);
k_lowSR_PAC2 = PL_lowSR.mu0*PL_lowSR.lambda^(PL_lowSR.n-1); % Ordinary PL coefficient
k_highSR_PAC2 = PL_highSR.mu0*PL_highSR.lambda^(PL_highSR.n-1); % Ordinary PL coefficient
n_PAC2=PL.n;
n_lowSR_PAC2 = PL_lowSR.n;
n_highSR_PAC2 = PL_highSR.n;

% Carreau-Yasuda
CarreauYasuda = createFit_CarreauYasuda(pac2_1(1:25),pac2_mean);
lambda_CarreauYasuda_PAC2 = CarreauYasuda.lambda;
close gcf;

%% Write model coefficients

% Model coefficient matrix for Table in paper
ModelCoefficients_PAC2 = { num2str(PL_highSR.mu0, 3) 'n/a' num2str(PL_highSR.lambda,3) num2str(PL_highSR.n,3) num2str(PL_highSR_gof.rsquare,3) num2str(PL_highSR_gof.sse,3);...
     num2str(Carreau.mu_0,3) num2str(Carreau.mu_inf, 3) num2str(Carreau.lambda,3) num2str(Carreau.n,3) num2str(Carreau_gof.rsquare,3) num2str(Carreau_gof.sse,3);...
     num2str(Cross.mu_0, 3) num2str(Cross.mu_inf, 3) num2str(Cross.lambda,3) num2str(Cross.n,3) num2str(Cross_gof.rsquare,3) num2str(Cross_gof.sse,3)};


parentpath = cd(cd('..'));
path = [parentpath '\Timescales'];
save([path '\Timescales_FC_PAC2.mat'],...
    'lambda_PL_PAC2',...
    'lambda_PL_lowSR_PAC2',...
    'lambda_PL_highSR_PAC2',...
    'lambda_Cr_PAC2',...
    'lambda_Ca_PAC2');

path = [parentpath '\NormalStressDifferences'];
save([path '\FC_Cross_PAC2.mat'],...
    'mu_0_PAC2',...
    'mu_inf_PAC2',...
    'lambda_Cr_PAC2',...
    'n_Cr_PAC2');

path = [parentpath '\FS'];
save([path '\FC_Cross_PAC2.mat'],...
    'mu_0_PAC2',...
    'mu_inf_PAC2',...
    'lambda_Cr_PAC2',...
    'n_Cr_PAC2');

save([path '\FC_PL_PAC2.mat'],...
    'n_PAC2',...
    'n_lowSR_PAC2',...
    'n_highSR_PAC2',...
    'k_PAC2',...
    'k_lowSR_PAC2',...
    'k_highSR_PAC2');

path = [parentpath '\Timescales'];
save([path '\FC_Cross_PAC2.mat'],...
    'mu_0_PAC2',...
    'mu_inf_PAC2',...
    'lambda_Cr_PAC2',...
    'n_Cr_PAC2');



%% PAC 4 - Laboratory 2 - UiS, Stavanger
addpath([pwd '\Data\Lab2_UiS'])
load Lab2;
MS = 8; % MarkerSize
MFC = 'k'; % MarkerFaceColor
AppVis=Lab2(2:2:end,:);
AppVis=AppVis([4:5],:);
AppVis = mean(AppVis,1);
SR = Lab2(7,:);
Lab2_PAC4_M = plot(SR,AppVis,...
    'o','Color',col_PAC4,'markersize',MS,'MarkerFaceColor','w');

%% Fit models
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\viscosity-based'
addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\viscosity-based'

% Cross
[Cross Cross_gof] = createFit_Cross(SR,AppVis);
close gcf;
plot(SR,Cross.mu_inf+(Cross.mu_0-Cross.mu_inf)./(1+(Cross.lambda.*SR).^Cross.n));
mu_0_PAC4 = Cross.mu_0;
mu_inf_PAC4 = Cross.mu_inf;
lambda_Cr_PAC4 = Cross.lambda;
n_Cr_PAC4 = Cross.n;

% Carreau
[Carreau Carreau_gof] = createFit_Carreau(SR,AppVis);
close gcf;
plot(SR,Carreau.mu_inf+(Carreau.mu_0-Carreau.mu_inf).*(1+(Carreau.lambda.*SR).^2).^((Carreau.n-1)./2));
lambda_Ca_PAC4 = Carreau.lambda;

% PL - Time constant version
[PL PL_gof] = createFit_PowerLaw_TimeConstantVersion(SR,AppVis,Carreau.mu_0);
close gcf;
plot(SR,PL.mu0.*(PL.lambda.*SR).^(PL.n-1));

[PL_lowSR PL_lowSR_gof] = createFit_PowerLaw_TimeConstantVersion(SR([1:8,33:40]),AppVis([1:8,33:40]),Carreau.mu_0);
close gcf;
plot(SR,PL_lowSR.mu0.*(PL_lowSR.lambda.*SR).^(PL_lowSR.n-1));

[PL_highSR PL_highSR_gof] = createFit_PowerLaw_TimeConstantVersion(SR([15:20,21:26]),AppVis([15:20,21:26]),Carreau.mu_0);
close gcf;
plot(SR,PL_highSR.mu0.*(PL_highSR.lambda.*SR).^(PL_highSR.n-1));

% Time scales
lambda_PL_PAC4 = PL.lambda;
lambda_PL_lowSR_PAC4 = PL_lowSR.lambda;
lambda_PL_highSR_PAC4 = PL_highSR.lambda;
    
% Ordinary PL coefficients
k_PAC4 = PL.mu0*PL.lambda^(PL.n-1);
k_lowSR_PAC4 = PL_lowSR.mu0*PL_lowSR.lambda^(PL_lowSR.n-1); % Ordinary PL coefficient
k_highSR_PAC4 = PL_highSR.mu0*PL_highSR.lambda^(PL_highSR.n-1); % Ordinary PL coefficient
n_PAC4=PL.n;
n_lowSR_PAC4 = PL_lowSR.n;
n_highSR_PAC4 = PL_highSR.n;

% Carreau-Yasuda
CarreauYasuda = createFit_CarreauYasuda(SR,AppVis);
lambda_CarreauYasuda = CarreauYasuda.lambda;
close gcf;


%% Uncertainty

% % 95% confidence bounds correspond to roughly 2SE
% % https://stats.stackexchange.com/questions/56596/finding-uncertainty-in-coefficients-from-polyfit-in-matlab
% level = 2*tcdf(-1,Carreau_gof.dfe);
% ci1 = confint(Carreau,level);
% SE1 = coeff-ci1(1,:);
% 
% % Relation of SD and SE
% % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1255808/
% SD = sqrt(13)*SE1;
% 





%% Write model coefficients

% Model coefficient matrix for Table in paper
ModelCoefficients_PAC4 = { num2str(PL_highSR.mu0, 3) 'n/a' num2str(PL_highSR.lambda,3) num2str(PL_highSR.n,3) num2str(PL_highSR_gof.rsquare,3) num2str(PL_highSR_gof.sse,3);...
     num2str(Carreau.mu_0, 3) num2str(Carreau.mu_inf, 3) num2str(Carreau.lambda,3) num2str(Carreau.n,3) num2str(Carreau_gof.rsquare,3) num2str(Carreau_gof.sse,3);...
     num2str(Cross.mu_0, 3) num2str(Cross.mu_inf, 3) num2str(Cross.lambda,3) num2str(Cross.n,3) num2str(Cross_gof.rsquare,3) num2str(Cross_gof.sse,3)};
path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Resub1\Tables';
xlswrite([path '\ModelCoefficients'],[ModelCoefficients_PAC2; ModelCoefficients_PAC4]);

parentpath = cd(cd('..'));
path = [parentpath '\Timescales'];
save([path '\Timescales_FC_PAC4.mat'],...
    'lambda_PL_PAC4',...
    'lambda_PL_lowSR_PAC4',...
    'lambda_PL_highSR_PAC4',...
    'lambda_Cr_PAC4',...
    'lambda_Ca_PAC4');

path = [parentpath '\NormalStressDifferences'];
save([path '\FC_Cross_PAC4.mat'],...
    'mu_0_PAC4',...
    'mu_inf_PAC4',...
    'lambda_Cr_PAC4',...
    'n_Cr_PAC4');

path = [parentpath '\FS'];
save([path '\FC_Cross_PAC4.mat'],...
    'mu_0_PAC4',...
    'mu_inf_PAC4',...
    'lambda_Cr_PAC4',...
    'n_Cr_PAC4');

save([path '\FC_PL.mat'],...
    'n_PAC4',...
    'n_lowSR_PAC4',...
    'n_highSR_PAC4',...
    'k_PAC4',...
    'k_lowSR_PAC4',...
    'k_highSR_PAC4');

path = [parentpath '\Timescales'];
save([path '\FC_Cross_PAC4.mat'],...
    'mu_0_PAC4',...
    'mu_inf_PAC4',...
    'lambda_Cr_PAC4',...
    'n_Cr_PAC4');


