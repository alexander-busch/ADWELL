% Clean up
clear all;
close all;
clc;

addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\export_fig');
% addpath('O:\OneDrive_NTNU\SWwd\MATLAB\Generic\export_fig');

% todo = 1; % read EXP data and setup CFD cases
todo = 0; % read CFD data and plot all results

v_ini = 0; % Compute vx0 and vy0 from experimental data
% v_ini = 1; % As given in paper

% Path of Fluent exports folder
% 'E:\... = Laptop
% 'O:\... = Desktop
resultspath = 'E:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - 3D\archive\171006_case_3D_C_BulkVelProf_noslip_Cross_SR=1_VM_PressGrad\exports';
% resultspath = 'E:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - 2D\archive\170903_2D_C_expVelProf_slip_Cross_SR=1\exports';

CFDcase = 'A';


% Definitions & parameters
Definitions;
Parameters;

% Tables
if todo == 1
    Cases = table;
    EXP = table;
    CFD_Khatibi = table;
else
    CFD_asis = table;
    CFD_new = table;
end

% Current row
row = 1;

% Build tables, loop all fluids, particle diameters and bulk velocities
for i = 1:length(fluidlist)
    for j = 1:length(d_p)
        for k = 1:length(U)
            
            if todo == 1
            
                % Build Cases table
                Cases.fluid(row,1) = string(fluidlist{i});
                Cases.d_p(row,1) = d_p(j);
                Cases.U(row,1) = U(k);
            
                % Build EXP tables
                EXPdata2Table;
            
            else
                
                % Build CFD tables
                CFDdata2Table;
                
            end
            
            % Increase table row index
            row = row + 1;
        end % of k
    end % of j
end % of i


if todo == 1
    % Delete empty rows = non-existing experiments 
    toDelete = EXP.vy0 == 0;
    EXP(toDelete,:) = [];
    Cases(toDelete,:) = [];
    CFD_Khatibi(toDelete,:) = [];

    % Get velocity profiles and write to EXP table
    % Loop all fluids, bulk velocities and particle diameters
    i_max = size(Cases);
    type = {'EXP', 'CFD'};
    for a = 1:2
          
        % read velocity profiles
        data = xlsread('exp\VelocityProfiles_Khatibi_et_al_(2016).xlsx', type{a},'B8:G283');
          
        for i = 1:i_max(1)
            if Cases.fluid(i,1) == 'H2O'
                if Cases.U(i) == 0.048
                    % Do nothing, no velocity profiles available
                else
                    if a == 1
                        EXP.y(i,1) = {data(:,1)};
                        EXP.u(i,1) = {data(:,2)};
                    else
                        CFD_Khatibi.y(i,1) = {data(:,1)};
                        CFD_Khatibi.u(i,1) = {data(:,2)};
                    end
                end    
            elseif Cases.fluid(i,1) == 'PAC2'
                if Cases.U(i) == 0.048
                    if a == 1
                        EXP.y(i,1) = {data(:,5)};
                        EXP.u(i,1) = {data(:,6)};
                    else
                        CFD_Khatibi.y(i,1) = {data(:,5)};
                        CFD_Khatibi.u(i,1) = {data(:,6)};
                    end
                else
                    EXP.y(i,1) = {0};
                    EXP.u(i,1) = {0};
                end
            elseif Cases.fluid(i,1) == 'PAC4'
                if Cases.U(i) == 0.048
                    EXP.y(i,1) = {0};
                    EXP.u(i,1) = {0};
                else
                    if a == 1
                        EXP.y(i,1) = {data(:,3)};
                        EXP.u(i,1) = {data(:,4)};
                    else
                        CFD_Khatibi.y(i,1) = {data(:,3)};
                        CFD_Khatibi.u(i,1) = {data(:,4)};
                    end
                end
            else
            end

        end
    end

    % Fit experimental velocity profiles and get slip velocities
    x0 = 0.02;
    y0 = 0.00;
    fig_VelocityProfiles = CreateFigure( 'Velocity profiles', 'u_x(y) [m/s]', 'y [m]' );
    for i=1:i_max(1)
        
        if EXP.y{i,1} ~= 0

            x = EXP.y{i,1} - x0;
            y = EXP.u{i,1};

            if Cases.fluid(i,1) == 'H2O'
                Fit = CreateVelocityProfileFit_H2O(x,y);
                EXP.u_slip(i,1) = Fit.p6;
            else
                Fit = CreateVelocityProfileFit_PAC(x,y);
                EXP.u_slip(i,1) = Fit.p3;
            end

            x_plot = linspace(-x0,0);

            if Cases.fluid(i,1) == 'H2O'
                plot(Fit.p1.*x_plot.^5 + Fit.p2.*x_plot.^4 + Fit.p3.*x_plot.^3 + Fit.p4.*x_plot.^2 + Fit.p5.*x_plot + Fit.p6,x_plot+x0,'Color',colorlist{1},'LineStyle','-');
                scatter(y,x+x0,'MarkerEdgeColor',colorlist{1},'Marker','o');
            else
                
                if Cases.fluid(i,1) == 'PAC2'
                    plot(Fit.p1.*x_plot.^2 + Fit.p2.*x_plot + Fit.p3,x_plot+x0,'Color',colorlist{2},'LineStyle','-');
                    scatter(y,x+x0,'MarkerEdgeColor',colorlist{2},'Marker','o');
                else
                    plot(Fit.p1.*x_plot.^2 + Fit.p2.*x_plot + Fit.p3,x_plot+x0,'Color',colorlist{3},'LineStyle','-');
                    scatter(y,x+x0,'MarkerEdgeColor',colorlist{3},'Marker','o');
                end               
                
            end

        end
    end
    EXP.u_slip(5,1) = Cases.U(5)/Cases.U(6)*EXP.u_slip(6,1);
    EXP.u_slip(7,1) = EXP.u_slip(5,1) ;

	% Save data
    save('Khatibi_et_al_2016.mat',...
    'Cases',...
    'EXP',...
    'CFD_Khatibi');
    

else
    % Delete empty rows = non-existing experiments 
    toDelete = CFD_new.delete == 1;
    CFD_asis(toDelete,:) = [];
    CFD_new(toDelete,:) = [];

    save('CFD.mat',...
        'CFD_asis',...
        'CFD_new');



    %% Plot all results

    
    load('Khatibi_et_al_2016.mat');
    load('CFD.mat');
    
    i_max = size(Cases);

    fig_titles = 'CFD vs. Exp. results';
    label_x = {'Fluid x-velocity u_x [m/s]'; 'x-coordinate [m]'};
    label_y = 'y-coordinate [m]';
    sub_fontsize = 10;
    
    fig = CreateFigure( fig_titles, '', '' );
    axis off;
    addpath E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files;
    fig_sub = tight_subplot(2,3,[.1 .06],[.06 .03],[.06 .02]);

    % Loop all fluids, bulk velocities and particle diameters and plot
    % results
    for i = 1:i_max(1) % i = 2, i = 3, i = 4
        
        % Subplot
        if Cases.fluid(i,1) == 'H2O'           
            subplot(fig_sub(2)); hold on;
            plotcolor = colorlist{1};
        elseif Cases.fluid(i,1) == 'PAC2'         
            subplot(fig_sub(3)); hold on;
            plotcolor = colorlist{2};
        elseif Cases.fluid(i,1) == 'PAC4'
            plotcolor = colorlist{3};
            if Cases.U(i) == 0.048
                subplot(fig_sub(5)); hold on;
            else
                subplot(fig_sub(6)); hold on;
            end
        else
        end
                    
        % Define markersize (= particle size)
        if Cases.d_p(i) == 0.00116
            markersize = markersizelist(1);
        elseif Cases.d_p(i) == 0.002
            markersize = markersizelist(2);
        else
            markersize = markersizelist(3);
        end

        % Check the length of vectors to be plotted as sometimes Fluent exports
        % inconsistent data in terms of array elements
        length_X = length(CFD_asis.X{i,1});
        length_Y = length(CFD_asis.Y{i,1});

        if length_X < length_Y
            CFD_asis.X{i,1}(length_X+1) = CFD_asis.X{i,1}(length_X);
            CFD_asis.t{i,1}(length_X+1) = CFD_asis.t{i,1}(length_X);    

        elseif length_X > length_Y
            CFD_asis.Y{i,1}(length_Y+1) = CFD_asis.Y{i,1}(length_Y);
        else
        end
            
        length_X = length(CFD_new.X{i,1});
        length_Y = length(CFD_new.Y{i,1});

        if length_X < length_Y
            CFD_new.X{i,1}(length_X+1) = CFD_new.X{i,1}(length_X);
            CFD_new.t{i,1}(length_X+1) = CFD_new.t{i,1}(length_X);    

        elseif length_X > length_Y
            CFD_new.Y{i,1}(length_Y+1) = CFD_new.Y{i,1}(length_Y);
        else
        end

        plot(CFD_asis.X{i,1}-CFD_asis.X{i,1}(1),CFD_asis.Y{i,1},'Color',plotcolor,'LineStyle','--'); % 1 is the row of the table and represents a particular case
        plot(CFD_new.X{i,1}-CFD_new.X{i,1}(1),CFD_new.Y{i,1},'Color',plotcolor,'LineStyle','-'); % 1 is the row of the table and represents a particular case
        plot(CFD_Khatibi.X{i,1},CFD_Khatibi.Y{i,1},'Color',plotcolor,'LineStyle',':'); % 1 is the row of the table and represents a particular case
        scatter(EXP.X{i,1},EXP.Y{i,1},markersize,'MarkerEdgeColor',plotcolor); % 1 is the row of the table and represents a particular case
        
        % Format subplot
        sub_title = [Cases.fluid{i,1} ', U = ' num2str(Cases.U(i,1)) 'm/s'];
        grid on;
        box on;
        title(sub_title,'FontSize',sub_fontsize);
        xlabel(label_x{2},'FontSize',sub_fontsize);
        ylabel(label_y,'FontSize',sub_fontsize);
        set(gca,...
            'FontSize',sub_fontsize,...
            'yLim',[0 0.02]);

    end % of i

    % Plot velocity profiles
    subplot(fig_sub(4)); hold on;
    rows = [1, 3, 6];
    legendentry = cell(length(rows)*4,1);
    for i=1:length(rows)
        
        % Experimental results and CFD results of Khatibi et al .(2016)
        scatter(EXP.u{rows(i),1},EXP.y{1,1},'MarkerEdgeColor',colorlist{i},'Marker','o');
        plot(CFD_Khatibi.u{rows(i),1},CFD_Khatibi.y{rows(i),1},'Color',colorlist{i},'LineStyle',':'); % 1 is the row of the table and represents a particular case
        
        % CFD results this study
        plot(CFD_asis.u_xy{rows(i),1},CFD_asis.y{rows(i),1},'Color',colorlist{i},'LineStyle','--'); % 1 is the row of the table and represents a particular case
        plot(CFD_new.u_xy{rows(i),1},CFD_new.y{rows(i),1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
        
        % plot(CFD_asis.u_xy_262mm{rows(i),1},CFD_asis.y_262mm{rows(i),1},'Color',colorlist{i},'LineStyle','--'); % 1 is the row of the table and represents a particular case
        % plot(CFD_new.u_xy_262mm{rows(i),1},CFD_new.y_262mm{rows(i),1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
        % 
        % plot(CFD_asis.u_xy_285mm{rows(i),1},CFD_asis.y_285mm{rows(i),1},'Color',colorlist{i},'LineStyle','--'); % 1 is the row of the table and represents a particular case
        % plot(CFD_new.u_xy_285mm{rows(i),1},CFD_new.y_285mm{rows(i),1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
        
        legendentry{i*4-3} = [ Cases.fluid{rows(i),1} ', EXP, U = ' num2str(Cases.U(rows(i))) ' m/s'];
        legendentry{i*4-2} = [ Cases.fluid{rows(i),1} ', CFD_K_h_a_t_i_b_i___e_t_a_l_._(_2_0_1_6_), U = ' num2str(Cases.U(rows(i))) ' m/s'];
        legendentry{i*4-1} = [ Cases.fluid{rows(i),1} ', CFD_a_s_i_s, U = ' num2str(Cases.U(rows(i))) ' m/s'];
        legendentry{i*4-0} = [ Cases.fluid{rows(i),1} ', CFD_n_e_w, U = ' num2str(Cases.U(rows(i))) ' m/s'];

    end
    % Plot PAC4, U = 0.048 m/s
    plot(CFD_asis.u_xy{5,1},CFD_asis.y{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    plot(CFD_new.u_xy{5,1},CFD_new.y{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    % plot(CFD_asis.u_xy_262mm{5,1},CFD_asis.y_262mm{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    % plot(CFD_new.u_xy_262mm{5,1},CFD_new.y_262mm{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    % plot(CFD_asis.u_xy_285mm{5,1},CFD_asis.y_285mm{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    % plot(CFD_new.u_xy_285mm{5,1},CFD_new.y_285mm{5,1},'Color',colorlist{i},'LineStyle','-'); % 1 is the row of the table and represents a particular case
    
    % Format subplot
    sub_title = 'Fluid velocity profiles';
    grid on;
    box on;
    title(sub_title,'FontSize',sub_fontsize);
    xlabel(label_x{1},'FontSize',sub_fontsize);
    ylabel(label_y,'FontSize',sub_fontsize);
    set(gca,...
        'FontSize',sub_fontsize,...
        'xLim',[0 0.16],...
        'yLim',[0 0.02]);
    % legend(legendentry,'location','west');
    
    
    % Case information
    subplot(fig_sub(1));
    if strcmp(CFDcase, 'A')
        sub_title = 'CFD case C (3D)';
    else
        sub_title = 'CFD case C (2D)';
    end
    title(sub_title,'FontSize',sub_fontsize);
    % Mesh
    inc=0.5;
    xmin=0;
    ymin=4;
    xmax=11;
    ymax=8;
    xsize = xmin:inc:xmax;
    ysize = ymin:inc:ymax;
    for i = 1:length(xsize)
        line([xsize(i) xsize(i)],[ysize(1) ysize(end)],'color','k');
    end
    for i = 1:length(ysize)
        line([xsize(1) xsize(end)],[ysize(i) ysize(i)],'color','k');
    end
    
    % Traps
    if strcmp(CFDcase, 'A')
        
        % First trap
        line([1 1],[2 4],'color','k');
        line([1.5 1.5],[2 4],'color','k');
        line([2 2],[2 4],'color','k');
        line([1 2],[3.5 3.5],'color','k');
        line([1 2],[3 3],'color','k');
        line([1 2],[2.5 2.5],'color','k');
        line([1 2],[2 2],'color','k');
        
        % Second trap
        line([2.5 2.5],[2 4],'color','k');
        line([3 3],[2 4],'color','k');
        line([3.5 3.5],[2 4],'color','k');
        line([2.5 3.5],[3.5 3.5],'color','k');
        line([2.5 3.5],[3 3],'color','k');
        line([2.5 3.5],[2.5 2.5],'color','k');
        line([2.5 3.5],[2 2],'color','k');
        
        % Third trap
        line([4 4],[2 4],'color','k');
        line([4.5 4.5],[2 4],'color','k');
        line([5 5],[2 4],'color','k');
        line([4 5],[3.5 3.5],'color','k');
        line([4 5],[3 3],'color','k');
        line([4 5],[2.5 2.5],'color','k');
        line([4 5],[2 2],'color','k');
    end
    
    set(gca,...
        'FontSize',sub_fontsize,...
        'xLim',[0 xmax+1],...
        'yLim',[0 xmax+1]);
    axis off;

    % x-axis
    text(xmax+1,ymin,'x','HorizontalAlignment' ,'left','VerticalAlignment','Top');
    arrow([xmax ymin],[xmax+1  ymin]);
    
    % y-axis
    text(0,ymax+2,'y','HorizontalAlignment' ,'right','VerticalAlignment','Bottom');
    arrow([0 ymax],[0 ymax+2]);
      
    % z-axis
    if strcmp(CFDcase, 'A')
        arrow([0 4],[-1 3]);
        text(-1,3,'z','HorizontalAlignment','right','VerticalAlignment','Middle');
    end
    
    % Particle injection
    text(1,10,'v_0 acc. to Khatibi et al. (2016)','Color','k','HorizontalAlignment' ,'left','VerticalAlignment','Bottom');
    arrow([1 10],[1 8]);
    
    % Velocity profiles
    text(xmin-0.2,0.5*(ymin+ymax),{'Bulk velocity acc. to', 'Khatibi et al. (2016)'},'Rotation',90,'HorizontalAlignment' ,'center','VerticalAlignment','Bottom');
    
    % Slip velocities
    if strcmp(CFDcase, 'A')
        text(0.5*(xmin+xmax),ymax,{'\color{black}u_s = 0 m/s'},'HorizontalAlignment' ,'left','VerticalAlignment','Bottom');
        text(0.5*(xmin+xmax),ymin,{'\color{black}u_s = 0 m/s'},'HorizontalAlignment' ,'left','VerticalAlignment','Top');
    else
        text(0.5*(xmin+xmax),ymax,{'\color{blue}u_s = 0.0218 m/s, \color{green}u_s = 0.002 m/s,','\color{red}u_s = 0.005 m/s, u_s = 0.009 m/s'},'HorizontalAlignment' ,'center','VerticalAlignment','Bottom');
        text(0.5*(xmin+xmax),ymin,{'\color{blue}u_s = 0.0240 m/s, \color{green}u_s = 0.011 m/s,','\color{red}u_s = 0.027 m/s, u_s = 0.050 m/s'},'HorizontalAlignment' ,'center','VerticalAlignment','Top');
    end
    
    % Shear rate model
    text(0,1,'Shear rate:','FontSize',9);
    text(2.5,1,{'$\dot{\gamma} = \frac{ || \mathbf{v_r} || }{d_p}$'},'interpreter','latex','FontSize',10,'HorizontalAlignment','Left','VerticalAlignment','Middle');
    % \sqrt{6}
    
    % Rheology
    % text(0,1,'Rheology: Carreau (1968), coefficients prev. slide','FontSize',9);
    text(0,0,'Rheology: Cross (1965), coefficients Figure 3. ','FontSize',9);
    
    % Drag model
    % text(0,0,'Drag: Chhabra & Uhlherr (1980)','FontSize',9');
    % text(0,0,'Drag: Acharya (1976)','FontSize',9');
    text(0,-1,'Drag: Schiller & Naumann (1933), Virtual mass and pressure gradient force','FontSize',9');
    
    % Save figures
    savefig(fig, [resultspath(1:length(resultspath)-7) 'ParticleTrajectories.fig']);
    set(gcf, 'color', 'none');
    set(findall(gcf,'type','axes'), 'color', 'none');
    export_fig x=f(y).png;
    movefile('x=f(y).png', resultspath(1:length(resultspath)-7));
	set(gcf, 'color', 'white');
	set(gca, 'color', 'white');

    % Save new CFD data
    save('Busch_et_al_2017.mat',...
        'CFD_asis',...
        'CFD_new');
    
    % Visualize velocity components of results
    velocity_components;
 
end