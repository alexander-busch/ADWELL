close all;
load('Khatibi_et_al_2016.mat');
    
rows = [1, 3, 6];

x0 = 0.02;
y0 = 0.00;

fig_VelocityProfiles = CreateFigure( 'Velocity profiles', 'u_x(y) [m/s]', 'y [m]' );
% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0.015;
factor_ver = 0.015;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];


for i=1:length(rows)

    x = EXP.y{rows(i),1} - x0;
    y = EXP.u{rows(i),1};
        
    if i == 1
        Fit = CreateVelocityProfileFit_H2O(x,y);
    else
        Fit = CreateVelocityProfileFit_PAC(x,y);
    end
    
    scatter(y,x+x0,'MarkerEdgeColor',colorlist{i},'Marker','o');
    x = linspace(-x0,0);

    if i == 1
        plot(Fit.p1.*x.^5 + Fit.p2.*x.^4 + Fit.p3.*x.^3 + Fit.p4.*x.^2 + Fit.p5.*x + Fit.p6,x+x0,'Color',colorlist{i},'LineStyle','-');
        str = {'H2O, U = 0.085 m/s',...
            ['u_x(y) = ' num2str(Fit.p1,2) ' Y^5 ' num2str(Fit.p2,2) ' Y^4 ' num2str(Fit.p3,2) ' Y^3 ' num2str(Fit.p4,2) ' Y^2 ' num2str(Fit.p5,2) ' Y +' num2str(Fit.p6,2)]};
        text(0.065,0.0188,str,'Color','blue','HorizontalAlignment' ,'left','VerticalAlignment','Bottom');
    else
        plot(Fit.p1.*x.^2 + Fit.p2.*x + Fit.p3,x+x0,'Color',colorlist{i},'LineStyle','-');
        if i == 2
            str = {'PAC2, U = 0.048 m/s',...
                ['u_x(y) = ' num2str(Fit.p1,2) ' Y^2 ' num2str(Fit.p2,2) ' Y +' num2str(Fit.p3,2)]};
            text(0.06,0.0038,str,'Color','green','HorizontalAlignment' ,'right','VerticalAlignment','Bottom');
        else
            str = {'PAC4, U = 0.085 m/s',...
                ['u_x(y) = ' num2str(Fit.p1,2) ' Y^2 ' num2str(Fit.p2,2) ' Y +' num2str(Fit.p3,2)]};
            text(0.116,0.0152,str,'Color','red','HorizontalAlignment' ,'left','VerticalAlignment','Bottom');
        end
        

    end

    
end



%% Print as image
addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\export_fig');
figurepath = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2017-xx - An Eulerian-Lagrangian CFD study of a particle settling in an orthogonal GNF shear flow\Figure\';

export_fig SlipVelocities.png;
movefile('SlipVelocities.png', figurepath);

% Specify Figure Size and Page Size
set(fig_VelocityProfiles,'PaperPositionMode','auto') %set paper pos for printing
fig_pos = fig_VelocityProfiles.PaperPosition;
fig_VelocityProfiles.PaperSize = [fig_pos(3) fig_pos(4)];

% Save Figure to File Format
figurepath = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2017-xx - An Eulerian-Lagrangian CFD study of a particle settling in an orthogonal GNF shear flow\Figure\';
fig_name = 'SlipVelocities';
print(fig_VelocityProfiles,[figurepath fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_VelocityProfiles,[figurepath fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_VelocityProfiles,[figurepath fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% saveas(fig, [path fig_name],'png') % save figure
% save2pdf([path fig_nam],fig,600) % requires m-file
