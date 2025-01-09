% Compute settling velocity
% accounting for particle-induced and fluctuation-induced shear and 
% corresponding effect on apparent viscosity.

%% Clean-up
clear all;
close all;
clc;
addpath 'E:\OneDrive - NTNU\SWwd\MATLAB\Generic\m-files';


%% Setup

% Parameters
Definitions;
Parameters;
ShearRateEstimate; close all; 
SR_mean = 0;

%% Settling velocity

% Create figures
CreateFigure('Drag coefficient', 'Re_p_0 [-]','c_D [-]');

% Create legend cell array
Legend = cell(length(fluidlist)+1,1);

% Plot Schiller-Naumann drag coefficient
Re_p = logspace(-5,4);
c_D = cD_SchillerNaumann1935( Re_p);
plot(Re_p,c_D,'--k');
Legend{1} = 'Schiller-Naumann (1935)';

% Plot Morrison drag coefficient
% c_D = cD_Morrison2016( Re_p);
% plot(Re_p,c_D,'--k');
% Legend{j+1} = 'Morrison ()';

% Create "Particle Reynolds numbers for zero flow" matrix
Re_p_0 = zeros(length(fluidlist),length(d_p));

% Loop shear rate range and compute c_D = f(Re_p) and v_set
for i = 1:length(SR_mean)
    
    % Loop all fluids
    for j = 1:length(fluidlist)
        
        % Fluid apparent viscosity based on wall shear rate estimate
        eta = CrossCoefficients.mu_inf{j}+(CrossCoefficients.mu_0{j}-CrossCoefficients.mu_inf{j})./(1+(CrossCoefficients.lambda{j}.*SR_mean(i)).^CrossCoefficients.n{j});
%         eta = CarreauCoefficients.mu_inf{j}+(CarreauCoefficients.mu_0{j}-CarreauCoefficients.mu_inf{j}).*((1+(CarreauCoefficients.lambda{j}.*SR_mean(i)).^2).^((CarreauCoefficients.n{j}-1)./2));

        
        
        
        % Fluctuating contribution to total shear rate
        SR_fluc = 0;
        SR_fluc = 0.01107539166.*SR_mean(i).*Re(i).^(5/16); % Acc. to Fluent R17.2 User Guide, page 287ff
        SR_fluc = 0.08333333333.*SR_mean(i).*Re(i).^(1/2); % Acc. to CFD online, https://www.cfd-online.com/Wiki/Introduction_to_turbulence/Turbulence_kinetic_energy
        
        % Iteratively compute settling velocity
        [ Re_p, c_D, v_set ] = IterateSettlingVelocity( d_p, rho_p, rho_f,...
            CrossCoefficients.mu_inf{j},...
            CrossCoefficients.mu_0{j},...
            CrossCoefficients.lambda{j},...
            CrossCoefficients.n{j},...
            SR_mean(i),...
            SR_fluc);
        
        
%          [ Re_p, c_D, v_set ] = IterateSettlingVelocity( d_p, rho_p, rho_f,...
%             CarreauCoefficients.mu_inf{j},...
%             CarreauCoefficients.mu_0{j},...
%             CarreauCoefficients.lambda{j},...
%             CarreauCoefficients.n{j},...
%             SR_mean(i),...
%             SR_fluc);
      
        % Save Particle Reynolds numbers for zero flow
        if i == 1
            Re_p_0(j,:) = Re_p;
        end % of if
        
        % Plot c_D vs. Re_p based on particle-induced shear rate
        % plot(Re_p,c_D,'o','MarkerEdgeColor',colorlist{j});
        scatter(Re_p_0(j,:),c_D,[30,80,200],'MarkerEdgeColor',colorlist{j}); % Re_p  
        
        % Create legend element
        Legend{j+1} = fluidlist{j};

    end % of for

end % of for

% Display legend
Legend = legend(Legend);
set(Legend,'Location','northeast');

% Clear workspace
clear Legend;


% Additional Formatting

% % Rectangle
% HSR = area([100 150], [ylim_max ylim_max]);
% HSR.FaceColor = [0.8 0.8 0.8];
% HSR.EdgeColor = 'none';
% HSR.FaceAlpha = 0.5;

% Text
% text(Re_p_0(5,1),c_D(1),['d_p = ' num2str(d_p(1)*1000) ' mm'],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top',...
%     'FontSize',12,'Rotation',-20); %'interpreter','latex'
% text(Re_p_0(5,2),c_D(2),['d_p = ' num2str(d_p(2)*1000) ' mm'],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top',...
%     'FontSize',12,'Rotation',-15); %'interpreter','latex'
% text(Re_p_0(5,3),c_D(3),['d_p = ' num2str(d_p(3)*1000) ' mm'],...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top',...
%     'FontSize',12,'Rotation',-10); %'interpreter','latex'
% 
% % Arrow
% annotation('textarrow',[0.225 0.225],[0.8 0.5],'String','SR')