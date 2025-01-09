% Plot flowcurves, rheometric data and Cross fit

%% Clean-up
clear all;
close all;
clc;
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files';
addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\viscosity-based';

%% Setup
Definitions;
CreateFigure('Flowcurves', 'Shear Rate [1/s]','Apparent viscosity [Pa.s]');

% Create legend cell array
Legend = cell(3*length(fluidlist),1);

% Create Shear rate vector for plotting Cross fit
SR = logspace(-2,4);

% Create table for model coefficients
Cross = cell2table(fluidlist','VariableNames',{'Fluid'});
Carreau = cell2table(fluidlist','VariableNames',{'Fluid'});
RheoData = cell2table(fluidlist','VariableNames',{'Fluid'});


%% Loop all fluids, get rheometric data, fit model, plot
for i = 1:length(fluidlist)
    
    % Read rheometric data (Shear rate & apparent viscosity) from Excel file
    if i ~= 3 % Not PAC2
        if i == 2
            Data = xlsread('Flowcurves.xlsx', fluidlist{i},'A7:B26');
            [num,source,raw] = xlsread('Flowcurves.xlsx', fluidlist{i},'A2');
            source = source{1};
        else
            Data = xlsread('Flowcurves.xlsx', fluidlist{i},'A5:B26');
            [num,source,raw] = xlsread('Flowcurves.xlsx', fluidlist{i},'A2');
            source = source{1};
        end
        
    else
        Data = [xlsread('Flowcurves.xlsx', fluidlist{i},'A5:B26'); xlsread('Flowcurves.xlsx', fluidlist{i},'D5:E26')];
        [num,source,raw] = xlsread('Flowcurves.xlsx', fluidlist{i},'A2');
        [num,source2,raw] = xlsread('Flowcurves.xlsx', fluidlist{i},'D2');        
        source = [source{1} ', ' source2{1}];
        clear source2;
    end % of if
    
    % Plot rheometric data
    scatter(Data(:,1), Data(:,2),'o','MarkerEdgeColor',colorlist{i});

    % Fit Cross model to rheometric data
    if i == 1 % H2O
        Cr.lambda = 0;
        Cr.mu_0 = 0.00102;
        Cr.mu_inf = 0.00102;
        Cr.n = 1;
    else % PAC
        Cr = createFit_Cross(Data(:,1),Data(:,2));
        close gcf;
    end % of if

    % Fit Carreau model to rheometric data
    if i == 1 % H2O
        Ca.lambda = 0;
        Ca.mu_0 = 0.00102;
        Ca.mu_inf = 0.00102;
        Ca.n = 1;
    else % PAC
        Ca = createFit_Carreau(Data(:,1),Data(:,2));
        close gcf;
    end % of if    
    
    % Plot Cross fit
    plot(SR,Cr.mu_inf+(Cr.mu_0-Cr.mu_inf)./(1+(Cr.lambda.*SR).^Cr.n),'--','Color',colorlist{i});
    
    % Plot Carreau fit
    plot(SR,Ca.mu_inf+(Ca.mu_0-Ca.mu_inf).*((1+(Ca.lambda.*SR).^2).^((Ca.n-1)./2)),'-','Color',colorlist{i});

    % Compute PL coefficents from Cross-fit
%     [n_PL, K_PL] = Cross2PL( Cross.lambda, Cross.n, Cross.mu_0, Cross.mu_inf, SR );
%     eta = K_PL.*SR.^(n_PL-1);
%     plot(SR, eta, 'b');
    
    % Write  model coefficients to table
    Cross.lambda{i} = Cr.lambda;
    Cross.mu_0{i} = Cr.mu_0;
    Cross.mu_inf{i} = Cr.mu_inf;
    Cross.n{i}= Cr.n;
    Cross.Source{i}= source;
    
    Carreau.lambda{i} = Ca.lambda;
    Carreau.mu_0{i} = Ca.mu_0;
    Carreau.mu_inf{i} = Ca.mu_inf;
    Carreau.n{i}= Ca.n;
    Carreau.Source{i}= source;
    
    RheoData.gamma_dot{i} = Data(:,1);
    RheoData.app_vis{i} = Data(:,2);
    
    % Create legend element
    Legend{3*i-2} = [fluidlist{i} '; Rheo. data of ' source];
    Legend{3*i-1} = [fluidlist{i} '; Cross fit (\lambda = ',...
        num2str(Cr.lambda,3) ', n = ',...
        num2str(Cr.n,3) ', \mu_0 = ',...
        num2str(Cr.mu_0,3) ', \mu_i_n_f = ',...
        num2str(Cr.mu_inf,3) ')'];
    Legend{3*i} = [fluidlist{i} '; Carreau fit (\lambda = ',...
        num2str(Ca.lambda,3) ', n = ',...
        num2str(Ca.n,3) ', \mu_0 = ',...
        num2str(Ca.mu_0,3) ', \mu_i_n_f = ',...
        num2str(Ca.mu_inf,3) ')'];
end

% Display legend
Legend = legend(Legend);
set(Legend,'Location','southwest');

% Plot PL fit of Khatibi et al. (2016)
%    plot(SR,0.082.*SR.^(0.78-1),'Color',colorlist{i});
%   plot(SR,0.256.*SR.^(0.71-1),'Color',colorlist{i});




%% Save results

% Model coefficients
save ('Cross.mat', 'Cross');
save ('Carreau.mat', 'Carreau');
save ('RheoData.mat', 'RheoData');

% Figure
savefig('Flowcurves.fig');
% close gcf;
% openfig('Flowcurves.fig');

