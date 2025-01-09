% clc;
% close all;
% clear all;

% Experimental and CFD data
% load('Khatibi_et_al_2016.mat');
% load('Busch_et_al_2017.mat');


% Length of caseslist
% i_max = size(Cases);

% Legend array
% legendentry = cell(2*i_max(1),1);

% Definitions & parameters
% Definitions;
% Parameters;


vel_comp = {'v_x'; 'v_y'};

for aa = 1:length(vel_comp)
    

    % Define figures
    fig_titles = {[vel_comp{aa}(3) ' = f(t)'];...
        ['d' vel_comp{aa}(3) '/dt = f(t)'];...
        ['d' vel_comp{aa}(3) '/dt = f(y)']};
    if strcmp(vel_comp{aa}, 'v_x')
        fig_titles{4} = ['d' vel_comp{aa}(3) '/dt / u(y) = f(y)'];
    end

    fig_axis = {[vel_comp{aa}(3) '-coordinate [m]'], 'Time [s]';...
        ['\Delta' vel_comp{aa}(3) '/\Deltat [m/s]'], 'Time [s]';...
        ['\Delta' vel_comp{aa}(3) '/\Deltat [m/s]'], 'y-coordinate [m]'};
    if strcmp(vel_comp{aa}, 'v_x') % Add normalized velocity figure
        fig_axis{4,1} = ['\Delta' vel_comp{aa}(3) '/\Deltat / u_f(y) [-]'];
        fig_axis{4,2} = 'y-coordinate [m]';
    else % Flip axis description
        fig_axis = fliplr(fig_axis);
    end
    sub_fontsize = 10;

    % Initialize figures and subplot matrices
    fig=zeros(length(fig_titles),1);
    fig_sub=zeros(length(fig),i_max(1));

    % Create index for looping subplots by rows first and then columns second
    index = reshape(1:8, 4, 2).';

    addpath E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files;

    % Create and format figures with subplots
    for i = 1:length(fig_titles) % Loop all figures

        fig(i) = CreateFigure( fig_titles{i}, '', '' );
        axis off;
        fig_sub(i,:) = tight_subplot(2,4,[.1 .06],[.06 .03],[.06 .02]);
        
        % Loop all cases and format subplots
        for j = 1:i_max(1) 

            % sub_title = [Cases.fluid{j,1} ', U = ' num2str(Cases.U(j,1)) 'm/s, d_p = ' num2str(Cases.d_p(j,1)*1000) ' mm'];
            % sub_fontsize = 8;
            % CreateSubplot( [2,4,index(j)], sub_title, fig_axis{i,1}, fig_axis{i,2}, sub_fontsize );
            
            subplot(fig_sub(i,index(j))); hold on;
            
            if i == 1    
                if strcmp(vel_comp{aa}, 'v_x')
                else
                    set(gca,'yLim',[0 0.02]);
                end   
            elseif i == 2
                if strcmp(vel_comp{aa}, 'v_x')

                else

                end  
            elseif i == 3
                if strcmp(vel_comp{aa}, 'v_x')
                    set(gca,'yLim',[0 0.02]);
                    
                    if j == 5 | j == 7
                        set(gca,'xLim',[0 0.15]);
                    elseif j == 6
                        set(gca,'xLim',[0 0.2]);
                    else
                    end
                    
                    
                    
                else
                    if j == 5 | j == 6
                        set(gca,...
                            'xLim',[0 0.02],...
                            'yLim',[-0.04 -0.01]);
                    elseif j == 7 | j == 8
                        set(gca,...
                            'xLim',[0 0.02],...
                            'yLim',[-0.07 -0.01]);
                    else
                        set(gca,...
                            'xLim',[0 0.02]);
                    end
                end  
            else
                set(gca,...
                    'xLim',[0 2],...
                    'yLim',[0 0.02]);
            end

            % Format subplot
            sub_title = [Cases.fluid{j,1} ', U = ' num2str(Cases.U(j,1)) 'm/s, d_p = ' num2str(Cases.d_p(j,1)*1000) ' mm'];
            grid on;
            box on;
            title(sub_title,'FontSize',sub_fontsize);
            xlabel(fig_axis{i,1},'FontSize',sub_fontsize);
            ylabel(fig_axis{i,2},'FontSize',sub_fontsize);
            set(gca,'FontSize',sub_fontsize);
            
        end

    end


    % Plot experimental data
    for j = 1:i_max(1)  

        % Fluid --> Color of plot symbols
        if Cases.fluid(j,1) == fluidlist{1} % 'H2O'
            plotcolor = colorlist{1};
        elseif Cases.fluid(j,1) == fluidlist{2} % 'PAC2'
            plotcolor = colorlist{2};
        elseif Cases.fluid(j,1) == fluidlist{3} % 'PAC4'
            plotcolor = colorlist{3}; 
        else
        end

        % Particle diameter --> Size of plot symbols
    %     if Cases.d_p(j) == d_p(1) % 0.00116
    %         markersize = markersizelist(1);
    %     elseif Cases.d_p(j) == d_p(2) % 0.002
    %         markersize = markersizelist(2);
    %     else
    %         markersize = markersizelist(3);
    %     end
        markersize = 5;

        % Get data from table
        t = EXP.t{j};
        x = EXP.X{j};
        y = EXP.Y{j};

        % Correct y-velocity vector for offset
        if y(1) > 0.02
            y = y-(y(1)-0.02);
        end

        % Smooth data in time (Average with span 5)
%         x_smooth = smooth(x);
%         y_smooth = smooth(y);
        
        x_smooth = smooth(x,9);
        y_smooth = smooth(y,9);
        
        x_smooth = smooth(x,9,'rloess');
        y_smooth = smooth(y,9,'rloess');
                
        % Horizontal velocity component
        dx = diff(x);
        dt = diff(t);
        dxdt = dx./dt;
        dx_smooth = diff(x_smooth);
        dxdt_smooth = dx_smooth./dt;

        % Vertical velocity component
        dy = diff(y);
        dt = diff(t);
        dydt = dy./dt;
        dy_smooth = diff(y_smooth);
        dydt_smooth = dy_smooth./dt;

        % diff() creates vectors with length n-1, so in order to plot
        % means of two consecutive entries are required  
        t_diff = 0.5 * (t(1:end-1) + t(2:end));
        x_diff = 0.5 * (x(1:end-1) + x(2:end));
        y_diff = 0.5 * (y(1:end-1) + y(2:end));
        x_smooth_diff = 0.5 * (x_smooth(1:end-1) + x_smooth(2:end));
        y_smooth_diff = 0.5 * (y_smooth(1:end-1) + y_smooth(2:end));
        
        for i = 1:length(fig_titles)

            figure(fig(i));
            subplot(fig_sub(i,index(j)));

            if i == 1
                if strcmp(vel_comp{aa}, 'v_x')
                    scatter(x,t,markersize,'MarkerEdgeColor',plotcolor);
                    set(gca,'YDir','reverse');
                else
                    scatter(t,y,markersize,'MarkerEdgeColor',plotcolor);
                    y_smooth = smooth(y);
                    
                end

                % legendentry{i} = ['Exp. (Khatibi et al. 2016), ', Cases.fluid{i,1}, ', U_f = ' num2str(Cases.U(i)) ' m/s, d_p = ' num2str(Cases.d_p(i)*1000) ' mm'];
            elseif i == 2
                if strcmp(vel_comp{aa}, 'v_x')
                    scatter(dxdt,t_diff,markersize,'MarkerEdgeColor',plotcolor);
                    plot(dxdt_smooth,t_diff,':','Color',plotcolor);
                    set(gca,'YDir','reverse');
                else
                    scatter(t_diff,dydt,markersize,'MarkerEdgeColor',plotcolor);
                    plot(t_diff,dydt_smooth,':','Color',plotcolor);
                end

            elseif i == 3
                if strcmp(vel_comp{aa}, 'v_x')
                    scatter(dxdt,y_diff,markersize,'MarkerEdgeColor',plotcolor);
                    plot(dxdt_smooth,y_smooth_diff,':','Color',plotcolor);
                    v_x_exp_smooth{j}=dxdt_smooth;
                    v_x_exp{j}=dxdt;
                else
                    scatter(y_diff,dydt,markersize,'MarkerEdgeColor',plotcolor);
                    plot(y_smooth_diff,dydt_smooth,':','Color',plotcolor);
                    set(gca,'XDir','reverse');
                    v_y_exp_smooth{j}=dydt_smooth;
                    v_y_exp{j}=dydt;
                    y_exp_smooth{j} = y_smooth_diff;
                    y_exp{j} = y_diff;
                end

            else

                % Fluid velocity profile u_f and y_f
                if j == 5 || j == 7 % PAC4, U_f = 0.048, no exp. velocity profiles availabe --> scale from 0.085
                    y_f = EXP.y{j+1,1};
                    u_f = Cases.U(j)./Cases.U(j+1).*EXP.u{j+1,1};
                else
                    y_f = EXP.y{j,1};
                    u_f = EXP.u{j,1};
                end

                % Remove NANs from fluid velocity profile u_f and y_f
                y_f(any(isnan(y_f),2),:)=[];
                u_f(any(isnan(u_f),2),:)=[];

                % Remove NANs from particle velocity dxdt and y
                % y(any(isnan(y),2),:)=[];
                y_diff(any(isnan(dxdt),2),:)=[];
                dxdt(any(isnan(dxdt),2),:)=[];

                % Fluid x-velocity at position of particle y
                u_f_y = interp1(y_f,u_f,y_diff,'spline');
                u_f_y_smooth = interp1(y_f,u_f,y_smooth_diff,'spline');
                
                scatter(dxdt./u_f_y,y_diff,markersize,'MarkerEdgeColor',plotcolor);
                plot(dxdt_smooth./u_f_y_smooth,y_smooth_diff,':','Color',plotcolor);
                
                u_x_exp_smooth{j}= u_f_y_smooth;
                u_x_exp{j}= u_f_y;
                
                
%                 figure; hold on;
%                 plot(u_f,y_f);
%                 plot(u_f_y,y_diff);
%                 plot(u_f_y_smooth,y_smooth_diff);
%                 
                
                
                
            end
        end 
    end


    % Plot CFD data
    for j = 1:i_max(1)  

        % Fluid --> Color of plot symbols
        if Cases.fluid(j,1) == fluidlist{1} % 'H2O'
            plotcolor = colorlist{1};
        elseif Cases.fluid(j,1) == fluidlist{2} % 'PAC2'
            plotcolor = colorlist{2};
        elseif Cases.fluid(j,1) == fluidlist{3} % 'PAC4'
            plotcolor = colorlist{3}; 
        else
        end

        % Get data from table
        t = CFD_new.t{j};
        x = CFD_new.X{j};
        y = CFD_new.Y{j};

        % Horizontal velocity component
        dx = diff(x);
        dt = diff(t);
        dxdt = dx./dt;

        % Vertical velocity component
        dy = diff(y);
        dt = diff(t);
        dydt = dy./dt;

        % diff() creates vectors with length n-1, so in order to plot
        % means of two consecutive entries are required  
        t_diff = 0.5 * (t(1:end-1) + t(2:end));
        x_diff = 0.5 * (x(1:end-1) + x(2:end));
        y_diff = 0.5 * (y(1:end-1) + y(2:end));

        for i = 1:length(fig_titles)

            figure(fig(i));
            subplot(fig_sub(i,index(j))); hold on;

            if i == 1
                if strcmp(vel_comp{aa}, 'v_x')
                    plot(x-x(1),t,'Color',plotcolor,'LineStyle','-');
                    % legendentry{i_max(1)+i} = ['CFD_n_e_w, ', Cases.fluid{i,1}, ', U_f = ' num2str(Cases.U(i)) ' m/s, d_p = ' num2str(Cases.d_p(i)*1000) ' mm'];
                else
                    plot(t,y,'Color',plotcolor,'LineStyle','-');
                end
            elseif i ==2
                if strcmp(vel_comp{aa}, 'v_x')
                    plot(dxdt,t_diff,'Color',plotcolor,'LineStyle','-');
                else
                    plot(t_diff,dydt,'Color',plotcolor,'LineStyle','-');
                end

            elseif i ==3
                if strcmp(vel_comp{aa}, 'v_x')
                    plot(dxdt,y_diff,'Color',plotcolor,'LineStyle','-');
                    v_x_cfd{j}=dxdt;
                else
                    plot(y_diff,dydt,'Color',plotcolor,'LineStyle','-');
                    v_y_cfd{j}=dydt;
                    y_cfd{j} = y_diff;
                end

            else
           
                
                
                
                % Fluid velocity profile u_f and y_f
                y_f = CFD_new.y{j,1};
                u_f = CFD_new.u_xy{j,1};

                % Remove NANs from fluid velocity profile u_f and y_f
                y_f(any(isnan(y_f),2),:)=[];
                u_f(any(isnan(u_f),2),:)=[];

                % Remove NANs from particle velocity dxdt and y
                % y(any(isnan(y),2),:)=[];
                y_diff(any(isnan(dxdt),2),:)=[];
                dxdt(any(isnan(dxdt),2),:)=[];

                % Fluid velocity at position of particle y
                u_f_y = interp1(y_f,u_f,y_diff,'spline');

                plot(dxdt./u_f_y,y_diff,'Color',plotcolor,'LineStyle','-');

                u_x_cfd{j}= u_f_y;
            end
        end 
    end
    
    
    
    % Schiller-Naumann settling velocities for steady-state particle settling
    if strcmp(vel_comp{aa}, 'v_y')
                
        load('Cross.mat')
        addpath('E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Particle settling\Particle settling in a shear flow');
        
        for j = 1:i_max(1)

            % Fluid --> Color of plot symbols
            if Cases.fluid(j,1) == 'H2O'
                plotcolor = colorlist{1};
            elseif Cases.fluid(j,1) == 'PAC2'
                plotcolor = colorlist{2};
            elseif Cases.fluid(j,1) == 'PAC4'
                plotcolor = colorlist{3}; 
            else
            end

            % Index of rheological coefficients table   
            if strcmp(Cases.fluid{j,1},'H2O')
                k=1;
            elseif strcmp(Cases.fluid{j,1},'PAC2')
                k=2;
            elseif strcmp(Cases.fluid{j,1},'PAC4')
                k=3;
            else
                
            end

            % Estimate Stokes settling velocity based on H2O [m/s]
            v_set = (rho_p-rho_f).*9.81.*(Cases.d_p(j)).^2./18./0.00102;
            v_set_old = 0;

            % Shear rate
            SR = v_set./(0.5.*Cases.d_p(j));

            % Accuracy
            accuracy = 1e-8;

            % Iteratively compute settling velocity
            while abs(v_set - v_set_old) > accuracy % eps(v_set)

                % Update settling velocity
                v_set_old = v_set;

                        % Compute settling velocity
                        v_set = SchillerNaumann1935(v_set,...
                            Cases.d_p(j),...
                            Cross.mu_inf{k},...
                            Cross.mu_0{k},...
                            Cross.lambda{k},...
                            Cross.n{k},...
                            rho_f,...
                            rho_p);

            end % of while
           
            for i = 2:length(fig_titles)

                figure(fig(i));
                subplot(fig_sub(i,index(j)));
            
                plot(xlim,[-v_set -v_set],'Color',plotcolor,'LineStyle','--');
                % legendentry{2*i_max(1)+i} = ['Schiller-Naumann, ', Cases.fluid{i,1}, ', U_f = ' num2str(Cases.U(i)) ' m/s, d_p = ' num2str(Cases.d_p(i)*1000) ' mm'];
            end
        end
    end
    
    % Swap position of subplot 6 and 7 to yield right format w.r.t. d_p and U_f
    for i = 1:length(fig_titles)
        pos1 = get(fig_sub(i,4) , 'position' );
        pos2 = get(fig_sub(i,7) , 'position' );
        set(fig_sub(i,4) , 'position' , pos2);
        set(fig_sub(i,7) , 'position' , pos1);
    end

        
	% Save figures
    savefig(fig, [resultspath(1:length(resultspath)-7) ['v_' vel_comp{aa} '.fig']]);
    set(fig, 'color', 'none');
    set(findall(fig,'type','axes'), 'color', 'none');
    for kk = 1:length(fig)
        figure(fig(kk));
        filename =  [vel_comp{aa} '_' num2str(kk) '.png'];
        expression = ['export_fig ' filename];
        eval(expression);
        newfilename = strrep(fig_titles{kk},'/','-');
        newfilename = strrep(newfilename,' ','');
        movefile(filename, [resultspath(1:length(resultspath)-7)  newfilename '.png']);
    end
    set(fig, 'color', 'white');
	set(findall(fig,'type','axes'), 'color', 'white');
end




fig = CreateFigure('Relative velocity ratios', '', '' );
axis off;
fig_sub=zeros(1,i_max(1));
fig_sub(1,:) = tight_subplot(2,4,[.1 .06],[.06 .03],[.06 .02]);
for j = 1:i_max(1)
        
    % Activate current supplot and plot velocity ratios    
    subplot(fig_sub(1,index(j))); hold on;
    
    % Format subplot
    sub_title = [Cases.fluid{j,1} ', U = ' num2str(Cases.U(j,1)) 'm/s, d_p = ' num2str(Cases.d_p(j,1)*1000) ' mm'];
    grid on;
    box on;
    title(sub_title,'FontSize',sub_fontsize);
    xlabel(fig_axis{3,1},'FontSize',sub_fontsize);
    ylabel('Rel. vel. ratio (v_x-u_x)/(v_y-u_y)','FontSize',sub_fontsize);
    set(gca,...
        'FontSize',sub_fontsize,...
        'xLim',[0 0.02]);
    if strcmp(Cases.fluid{j,1},'H2O')
        set(gca,'yLim',[-0.4 0.4]);
    else
        set(gca,'yLim',[-1 1]);
    end
        

    % Fluid --> Color of plot symbols
    if Cases.fluid(j,1) == fluidlist{1} % 'H2O'
        plotcolor = colorlist{1};
    elseif Cases.fluid(j,1) == fluidlist{2} % 'PAC2'
        plotcolor = colorlist{2};
    elseif Cases.fluid(j,1) == fluidlist{3} % 'PAC4'
        plotcolor = colorlist{3}; 
    else
    end
    
     % Remove NAN
    v_x_exp{j}(any(isnan(v_x_exp{j}),2),:)=[];
    v_y_exp{j}(any(isnan(v_y_exp{j}),2),:)=[];
    y_exp{j}(any(isnan(y_exp{j}),2),:)=[];
    y_exp{j}(y_exp{j}<0,:)=[];
    v_x_cfd{j}(any(isnan(v_x_cfd{j}),2),:)=[];
    v_y_cfd{j}(any(isnan(v_y_cfd{j}),2),:)=[];
    y_cfd{j}(any(isnan(y_cfd{j}),2),:)=[];
    
    % Correct length
    if length(y_exp{j})>length(v_y_exp{j})
        y_exp{j}(end)=[];
    elseif length(y_exp{j})<length(v_y_exp{j})
        v_x_exp{j}(end)=[];
        u_x_exp{j}(end)=[];
        v_y_exp{j}(end)=[];
    end
    
    % Plot relative velocity ratios                   
    plot(y_exp_smooth{j},(v_x_exp_smooth{j}-u_x_exp_smooth{j})./v_y_exp_smooth{j},'Color',plotcolor,'LineStyle',':');
    scatter(y_exp{j},(v_x_exp{j}-u_x_exp{j})./v_y_exp{j},markersize,'MarkerEdgeColor',plotcolor);
    plot(y_cfd{j},(v_x_cfd{j}-u_x_cfd{j})./v_y_cfd{j},'Color',plotcolor,'LineStyle','-');
    set(gca, 'Xdir', 'reverse');
        
    length(v_x_exp{j})
    length(u_x_exp{j})
    length(v_y_exp{j})
    length(y_exp{j})
end

% Swap position of subplot 6 and 7 to yield right format w.r.t. d_p and U_f
pos1 = get(fig_sub(1,4) , 'position' );
pos2 = get(fig_sub(1,7) , 'position' );
set(fig_sub(1,4) , 'position' , pos2);
set(fig_sub(1,7) , 'position' , pos1);



% Save figure
savefig(fig, [resultspath(1:length(resultspath)-7) 'RelVelRatio.fig']);
set(gcf, 'color', 'none');
set(findall(gcf,'type','axes'), 'color', 'none');
export_fig RelVelRatio.png;
movefile('RelVelRatio.png', resultspath(1:length(resultspath)-7));
set(gcf, 'color', 'white');
set(gca, 'color', 'white');

% figure(fig_x_t);
% legend(legendentry,'location','northeast');
% set(gca,'yLim',[0 0.02]);
% figure(fig_dxdt_t);
% % legend(legendentry,'location','southeast');
% figure(fig_dxdt_y);
% % legend(legendentry,'location','southeastoutside');
% set(gca,'xLim',[0 0.02]);


%% Print as image

% Expand Axes to Fill Figure (--> Minimum white space)
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% factor_hor = 0 %.003;
% factor_ver = 0 %.01;
% left = outerpos(1) + ti(1) + factor_hor;
% bottom = outerpos(2) + ti(2) + factor_ver;
% ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
% ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
% ax.Position = [left bottom ax_width ax_height];

% % Specify Figure Size and Page Size
% set(fig,'PaperPositionMode','auto') %set paper pos for printing
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% 
% % Save Figure to File Format
% path = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2016-05 - Rheological properties of PAC\Paper - V3 - Applied Rheology - Timescale Overview\Figures\';
% fig_name = 'FC';
% print(fig,[path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig,[path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(fig,[path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% % saveas(fig, [path fig_name],'png') % save figure
% % save2pdf([path fig_nam],fig,600) % requires m-file



% addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files');
% 
% jointfig(fig(i),2,4) 
% 
% 
% [ha, pos] = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01]) 
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end 
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')



%% Subplots and one legend

% Construct a figure with subplots and data
% fig = CreateFigure( 'y-component of particle velocity', '', '' );
% axis_des = {'Time [s]', '\Deltax/\Deltat [m/s]'};
% fig_sub_1 = CreateSubplot( [2,2,1], 'H2O, U_f = 0.048 m/s', axis_des{1}, axis_des{2} );
% fig_sub_2 = CreateSubplot( [2,2,2], 'PAC2, U_f = 0.048 m/s', axis_des{1}, axis_des{2} );
% fig_sub_3 = CreateSubplot( [2,2,3], 'PAC4, U_f = 0.048 m/s', axis_des{1}, axis_des{2} );
% fig_sub_4 = CreateSubplot( [2,2,4], 'PAC4, U_f = 0.085 m/s', axis_des{1}, axis_des{2} );
% 
% 
% 
% subplot(fig_sub_1)
% 
% subplot(fig_sub_2)
% 
% subplot(fig_sub_3)
% 
% subplot(fig_sub_4)







% 
% line1 = plot(1:10,rand(1,10),'b');
% title('Axes 1');
% subplot(2,2,2);
% line2 = plot(1:10,rand(1,10),'g');
% title('Axes 2');
% subplot(2,2,3);
% line3 = plot(1:10,rand(1,10),'r');
% title('Axes 3');
% a = subplot(2,2,4);
% line4 = plot(1:10,rand(1,10),'v_y');
% title('Axes 4');

% Construct a Legend with the data from the sub-plots
% hL = legend([line1,line2,line3,line4],{'Data Axes 1','Data Axes 2','Data Axes 3','Data Axes 4'});
% hL = legend([line1,line2,line3],{'Data Axes 1','Data Axes 2','Data Axes 3'});
% 
% test = get(a,'Position');
% delete(a)
% % Programatically move the Legend
% newPosition = test; % [0.4 0.4 0.3 0.3];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits);
% 
% 
% 
% fig_dxdt_y = CreateFigure( 'dx/dt = f(y)', 'y-coordinate [m]', '\Deltax/\Deltat [m/s]' );