% Create x-y-plots of bed height for various CFD cases

% Clean up
clear all;
close all;
clc;


% Path of paper figure folder
figurepath = 'M:\Documents\AdWell\8 - Publications & Conferences\2017-05 - A bed morphodynamics model\Paper\Figures';


% H2O - Usl = 0.26 - dp = 0.0003
figurename{1} = '\BedHeight_H2O_Usl0.26_dp0.0003';
path{1} = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170221_1 - MP + P + Debugging\test\170310 - H2O - Usl0.26 - dp0.0003\bed';
flowtimearray{1} = [0200, 0400, 0600, 0800, 1000, 2000, 3000, 4000, 5000, 6000, 7000];

% H2O - Usl = 0.43 - dp = 0.0003
figurename{2} = '\BedHeight_H2O_Usl0.43_dp0.0003';
path{2} = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170221_1 - MP + P + Debugging\test\170308 - H2O - Usl0.45 - d0.0003\bed';
flowtimearray{2} = [5000, 6000, 7000, 8000, 9000, 10000, 11000];

% H2O - Usl = 0.44 - dp = 0.0012
figurename{3} = '\BedHeight_H2O_Usl0.44_dp0.0012';
path{3} = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170221_1 - MP + P + Debugging\bed';
flowtimearray{3} = [0200, 0400, 0600, 0800, 1000];

% PAC1 - Usl = 0.45 - dp = 0.0012
figurename{4} = '\BedHeight_PAC1_Usl0.45_dp0.0012';
path{4} = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170310_1 - MP + P + nN\test\170310 - PAC1 - Usl0.45 - dp1.2\bed';
flowtimearray{4} = [0200, 0400, 0600, 0800, 1000, 2000, 3000, 4000, 5000, 6000];

% PAC1 - Usl = 0.81 - dp = 0.0012
figurename{5} = '\BedHeight_PAC1_Usl0.81_dp0.0012';
path{5} = 'C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170310_1 - MP + P + nN\test\170313 - PAC1 - Usl0.81 - dp1.2\bed';
flowtimearray{5} = [0200, 0400, 0600, 0800, 1000, 1100, 1200, 1300, 1400, 1500];



% Create all figures
for i = 1:length(path)
    % Create figure
    fig(i) = figure('color','w');
    hold on;
    grid on;
    N = length(flowtimearray{i}); % Color definition: 
    ColorSet = ColorBand(N);
    set(gca,...
        'xlim', [0 3],...
        'ylim', [0 0.04],...
        'box','on',...
        'ColorOrder', ColorSet);
    % 'FontSize',24)
    xlabel('x [m]');
    ylabel('y [m]');
    
    flowtime = flowtimearray{i};
    % Read and plot data
    for j = 1:N
        if flowtime(j)<1000
            filename=cat(2,'\movingbed_0',num2str(flowtime(j)),'_x_y.txt');
        else
            filename=cat(2,'\movingbed_',num2str(flowtime(j)),'_x_y.txt');
        end

        filepath = strcat(path(i),filename);
        txt=fileread(filepath{1});
        % Skip first column of numeric data, read two columns, then ignore the rest
        format=[repmat('%*f', [1 1]) '%f%f%*[^\n]'];  
        data = textscan(txt, format, 'Delimiter', '\t', 'HeaderLines', 1);
        x = data(:,1);
        y = data(:,2);
        x =x{1};
        y=y{1};
        plot(x,y,'linewidth',3);
        legend_str{j} = num2str(flowtime(j)/100','t = %-d s');
    end
    
   h_legend=columnlegend(N/2,legend_str,'Location','North');
    %set(h_legend,'FontSize',40);
   legend_str=cell(1);
    pathstring = strcat(figurepath, figurename(i));
    
    %set(findall(fig(i),'-property','FontSize'),'FontSize',24)
    
    
    % export_fig pathstring
    % print(fig, pathstring{1},'-dpng');  % '-dpng' '-depsc' 'jpeg' 'eps'   
    % saveas(gcf,pathstring{1},'pdf')
end
