close all;
clear all;
clc;

addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
addpath([pwd '\utilities']);
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Wellbore\Figures\';
load('wellboredata.mat');

%% Effect of plain rotation: CTR for H2O

% Select case2print, last gets indices of cell array elements w/ lengths > 1
cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'H2O'),:);

CreateFigure('CTR vs. rpm and \Deltap/\Deltax for H2O','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 0.045],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        Z(:,bb) = cases2print.CTR_Q{aa,bb}';
    end
    
    % Plot CTR
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),max(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),max(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('CTR [-]');

TightFigure(1);

% set annulus symbols
% xpos = .91;
% 
% for aa = 1:5
%     % aa = 4
%     if (cases2print.e(aa)<0 && aa<4)
%         imagename = 'e-.png';
%         imagepos = [xpos .20 .1 .1];
%     elseif (cases2print.e(aa)==0 && aa<4)
%         imagename = 'e0.png';
%         imagepos = [xpos .31 .1 .1];
%     elseif (cases2print.e(aa)>0 && aa<4)
%         imagename = 'e+.png';
%         imagepos = [xpos .43 .1 .1];
%     else
%         if aa==4
%             imagename = 'aw.png';
%             imagepos = [xpos .35 .1 .1];
%         else
%             imagename = 'sw.png';
%             imagepos = [xpos .55 .1 .1];
%         end
%         
%     end
% 
%   
%     if aa<4
%         [indexedImage, cmap]=imread(imagename);
%         trans = indexedImage~=0;
%         ax2(aa) = axes('pos',imagepos,'Color', 'none');
%         axes(ax2(aa));
%         f = imshow( indexedImage, cmap);
%         set(f, 'AlphaData', trans);
%     else
%         
%     end
% end

% Print to files
fig_name = 'MP_CTR-rpm-dpdx_H2O';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Effect of plain rotation: U_f for H2O

CreateFigure('U_f vs. rpm and \Deltap/\Deltax for H2O','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 2.5],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.U_f_x{aa,bb}>0
            Z(:,bb) = cases2print.U_f_x{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot U_f
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),min(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),min(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X))+5,max(max(Y))+10,max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('U_f_x [m/s]');

TightFigure(1);

% set annulus symbols
% xpos = .91;
% 
% for aa = 1:5
%     % aa = 4
%     if (cases2print.e(aa)<0 && aa<4)
%         imagename = 'e-.png';
%         imagepos = [xpos .20 .1 .1];
%     elseif (cases2print.e(aa)==0 && aa<4)
%         imagename = 'e0.png';
%         imagepos = [xpos .31 .1 .1];
%     elseif (cases2print.e(aa)>0 && aa<4)
%         imagename = 'e+.png';
%         imagepos = [xpos .43 .1 .1];
%     else
%         if aa==4
%             imagename = 'aw.png';
%             imagepos = [xpos .35 .1 .1];
%         else
%             imagename = 'sw.png';
%             imagepos = [xpos .55 .1 .1];
%         end
%         
%     end
% 
%   
%     if aa<4
%         [indexedImage, cmap]=imread(imagename);
%         trans = indexedImage~=0;
%         ax2(aa) = axes('pos',imagepos,'Color', 'none');
%         axes(ax2(aa));
%         f = imshow( indexedImage, cmap);
%         set(f, 'AlphaData', trans);
%     else
%         
%     end
% end

% Print to files
fig_name = 'MP_U_f-rpm-dpdx_H2O';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Effect of plain rotation: ROP for H2O

CreateFigure('ROP vs. rpm and \Deltap/\Deltax for H2O','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 120],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.ROP{aa,bb}>0
            Z(:,bb) = cases2print.ROP{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot ROP
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),max(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),max(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('ROP [m/h]');

TightFigure(1);

% Print to files
fig_name = 'MP_ROP-rpm-dpdx_H2O';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 


%% Effect of plain rotation: Re for H2O

CreateFigure('Re vs. rpm and \Deltap/\Deltax for H2O','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 200000],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.Re_MR{aa,bb}>0
            Z(:,bb) = cases2print.Re_MR{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot ROP
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),min(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),min(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),min(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('Re_M_R [-]');

TightFigure(1);

% Print to files
fig_name = 'MP_Re-rpm-dpdx_H2O';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 




%% Effect of plain rotation: CTR for PAC1

% Select case2print, last gets indices of cell array elements w/ lengths > 1
cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'PAC1'),:);

CreateFigure('U_f, U_s, vs. rpm and \Deltap/\Deltax for PAC1','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 0.045],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 2

    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        Z(:,bb) = cases2print.CTR_Q{aa,bb}';
    end
    
    % Plot CTR
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),max(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),max(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('CTR [-]');

TightFigure(1);

% set annulus symbols
% xpos = .91;
% 
% for aa = 1:5
%         % aa = 4
%     if (cases2print.e(aa)<0 && aa<4)
%         imagename = 'e-.png';
%         imagepos = [xpos .35 .1 .1];
%     elseif (cases2print.e(aa)==0 && aa<4)
%         imagename = 'e0.png';
%         imagepos = [xpos .49 .1 .1];
%     elseif (cases2print.e(aa)>0 && aa<4)
%         imagename = 'e+.png';
%         imagepos = [xpos .63 .1 .1];
%     else
%         if aa==4
%             imagename = 'aw.png';
%             imagepos = [xpos .35 .1 .1];
%         else
%             imagename = 'sw.png';
%             imagepos = [xpos .55 .1 .1];
%         end
%         
%     end
% 
%   
%     if aa<4
%         [indexedImage, cmap]=imread(imagename);
%         trans = indexedImage~=0;
%         ax2(aa) = axes('pos',imagepos,'Color', 'none');
%         axes(ax2(aa));
%         f = imshow( indexedImage, cmap);
%         set(f, 'AlphaData', trans);
%     else
%         
%     end
% end

% Print to files
fig_name = 'MP_CTR-rpm-dpdx_PAC1';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Effect of plain rotation: U_f for PAC1

CreateFigure('U_f vs. rpm and \Deltap/\Deltax for PAC1','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 2.5],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 4
    % aa = 5
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.U_f_x{aa,bb}>0
            Z(:,bb) = cases2print.U_f_x{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot U_f
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),min(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),min(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X))+2,max(max(Y))+20,max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('U_f_x [m/s]');

TightFigure(1);

% set annulus symbols
% xpos = .91;
% 
% for aa = 1:5
%     % aa = 4
%     if (cases2print.e(aa)<0 && aa<4)
%         imagename = 'e-.png';
%         imagepos = [xpos .20 .1 .1];
%     elseif (cases2print.e(aa)==0 && aa<4)
%         imagename = 'e0.png';
%         imagepos = [xpos .31 .1 .1];
%     elseif (cases2print.e(aa)>0 && aa<4)
%         imagename = 'e+.png';
%         imagepos = [xpos .43 .1 .1];
%     else
%         if aa==4
%             imagename = 'aw.png';
%             imagepos = [xpos .35 .1 .1];
%         else
%             imagename = 'sw.png';
%             imagepos = [xpos .55 .1 .1];
%         end
%         
%     end
% 
%   
%     if aa<4
%         [indexedImage, cmap]=imread(imagename);
%         trans = indexedImage~=0;
%         ax2(aa) = axes('pos',imagepos,'Color', 'none');
%         axes(ax2(aa));
%         f = imshow( indexedImage, cmap);
%         set(f, 'AlphaData', trans);
%     else
%         
%     end
% end

% Print to files
fig_name = 'MP_U_f-rpm-dpdx_PAC1';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'


%% Effect of plain rotation: ROP for PAC1

CreateFigure('ROP vs. rpm and \Deltap/\Deltax for PAC1','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 120],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.ROP{aa,bb}>0
            Z(:,bb) = cases2print.ROP{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot ROP
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),max(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),max(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),max(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),max(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('ROP [m/h]');

TightFigure(1);

% set annulus symbols
% xpos = .91;
% 
% for aa = 1:5
%     % aa = 4
%     if (cases2print.e(aa)<0 && aa<4)
%         imagename = 'e-.png';
%         imagepos = [xpos .20 .1 .1];
%     elseif (cases2print.e(aa)==0 && aa<4)
%         imagename = 'e0.png';
%         imagepos = [xpos .31 .1 .1];
%     elseif (cases2print.e(aa)>0 && aa<4)
%         imagename = 'e+.png';
%         imagepos = [xpos .43 .1 .1];
%     else
%         if aa==4
%             imagename = 'aw.png';
%             imagepos = [xpos .35 .1 .1];
%         else
%             imagename = 'sw.png';
%             imagepos = [xpos .55 .1 .1];
%         end
%         
%     end
% 
%   
%     if aa<4
%         [indexedImage, cmap]=imread(imagename);
%         trans = indexedImage~=0;
%         ax2(aa) = axes('pos',imagepos,'Color', 'none');
%         axes(ax2(aa));
%         f = imshow( indexedImage, cmap);
%         set(f, 'AlphaData', trans);
%     else
%         
%     end
% end

% Print to files
fig_name = 'MP_ROP-rpm-dpdx_PAC1';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 


%% Plot U_f and U_s vs. dpdx for 130 rpm and H2O, e+ vs sw

% Select case2print, last gets indices of cell array elements w/ lengths > 1
cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'H2O'),:);

fig = CreateFigure('U_f, U_s vs. dpdx for 130 rpm and H2O, e+ vs sw','','', {'lin' 'lin'}, 'DINA5');

% Select rpm
bb = 5; % 130 rpm

% Set axis limits
% set(gca,...
%     'xlim', [100 550],...
%     'ylim', [1e-3 2e0]);
set(gca,...
    'FontSize',12)

for aa = 1:5
    % aa = 2

    % Determine sign of plot eccentricity
    if (cases2print.e(aa)<0 && aa<4)
%         sign = -1;
        markersymbol = 'v';
    elseif (cases2print.e(aa)==0 && aa<4)
%         sign = 0;
        markersymbol = 'o';
    elseif (cases2print.e(aa)>0 && aa<4)
%         sign = 1;
        markersymbol = '^';
        plotcolor = 'k';
    else
        if aa==4
             markersymbol = '<';
        else
             markersymbol = '>';
             plotcolor = 'r';
        end
    end        
    
    if (aa == 3 || aa == 5)
        % Plot fluid superficial velocity vs dpdx
        %yyaxis left
        subplot(2,1,1); hold on;
        plot(-cases2print.dpdx{aa},cases2print.U_f_x{aa,bb},...
            'LineStyle','-',...
            'Color',plotcolor,...
            'Marker', '.',...
            'MarkerSize',10,...
            'MarkerEdgeColor',plotcolor,...
            'MarkerFaceColor',plotcolor);

        % Plot solid superficial velocity vs dpdx
        subplot(2,1,2); hold on;
        plot(-cases2print.dpdx{aa},cases2print.U_s_x{aa,bb},...
            'LineStyle','-',...
            'Color',plotcolor,...
            'Marker', '.',...
            'MarkerSize',10,...
            'MarkerEdgeColor',plotcolor,...
            'MarkerFaceColor',plotcolor);
    end

end

% set(gca,...
%     'yscale', 'log',...
%     'YColor','k');

% ylim
% ylabel('\fontsize{11}Time averaged {\color[rgb]{0 0.4470 0.7410}U_f}, {\color[rgb]{0.8500 0.3250 0.0980}U_s} [m/s]','interpreter','tex');
% 

for ii = 1:2
subplot(2,1,ii);
xlabel('Pressure gradient -\Deltap/\Deltax [Pa/m]');
if ii==1
    ylabel('(Time averaged) U_f [m/s]');
    set(gca,'Ylim',[0 1.5]);
else
    ylabel('(Time averaged) U_s [m/s]');
    set(gca,'Ylim',[0 0.05]);
end
box on;
grid on;
TightFigure(1);
end


% Print to files
fig_name = 'MP_U-dpdx_H2O';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Plot U_f and U_s vs. dpdx for 130 rpm and PAC1, e+ vs sw

% Select case2print, last gets indices of cell array elements w/ lengths > 1
cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'PAC1'),:);

fig = CreateFigure('U_f, U_s vs. dpdx for 130 rpm and PAC1, e+ vs sw','','', {'lin' 'lin'}, 'DINA5');

% Select rpm
bb = 5; % 130 rpm

% Set axis limits
% set(gca,...
%     'xlim', [100 550],...
%     'ylim', [1e-3 2e0]);

set(gca,...
    'FontSize',12)

for aa = 1:5
    % aa = 2

    % Determine sign of plot eccentricity
    if (cases2print.e(aa)<0 && aa<4)
%         sign = -1;
        markersymbol = 'v';
    elseif (cases2print.e(aa)==0 && aa<4)
%         sign = 0;
        markersymbol = 'o';
    elseif (cases2print.e(aa)>0 && aa<4)
%         sign = 1;
        markersymbol = '^';
        plotcolor = 'k';
    else
        if aa==4
             markersymbol = '<';
        else
             markersymbol = '>';
             plotcolor = 'r';
        end
    end        
    
    if (aa == 3 || aa == 5)
        % Plot fluid superficial velocity vs dpdx
        %yyaxis left
        subplot(2,1,1); hold on;
        
        
        U_r = 2.*pi.*cases2print.rpm{aa}./60.*cases2print.d_i(aa);
        
        plot(-cases2print.dpdx{aa},U_r(bb)./cases2print.U_f_x{aa,bb},...
            'LineStyle','-',...
            'Color',plotcolor,...
            'Marker', '.',...
            'MarkerSize',10,...
            'MarkerEdgeColor',plotcolor,...
            'MarkerFaceColor',plotcolor);

        % Plot solid superficial velocity vs dpdx
        subplot(2,1,2); hold on;
        plot(-cases2print.dpdx{aa},cases2print.U_s_x{aa,bb},...
            'LineStyle','-',...
            'Color',plotcolor,...
            'Marker', '.',...
            'MarkerSize',10,...
            'MarkerEdgeColor',plotcolor,...
            'MarkerFaceColor',plotcolor);
    end

end

% set(gca,...
%     'yscale', 'log',...
%     'YColor','k');

% ylim
% ylabel('\fontsize{11}Time averaged {\color[rgb]{0 0.4470 0.7410}U_f}, {\color[rgb]{0.8500 0.3250 0.0980}U_s} [m/s]','interpreter','tex');
% 

for ii = 1:2
subplot(2,1,ii);
xlabel('Pressure gradient -\Deltap/\Deltax [Pa/m]');
if ii==1
    ylabel('(Time averaged) U_f [m/s]');
    set(gca,'Ylim',[0 8]);
else
    ylabel('(Time averaged) U_s [m/s]');
    set(gca,'Ylim',[0 0.05]);
end
box on;
grid on;
TightFigure(1);
end


% Print to files
fig_name = 'MP_U-dpdx_PAC1';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Effect of plain rotation: Re for PAC1

CreateFigure('Re vs. rpm and \Deltap/\Deltax for PAC1','Rotation rate [rpm]','\Deltap/\Deltax  [Pa/m]', {'lin' 'lin'}, 'DINA5');
ax1 = gca;

% Set axis limits
set(ax1,...
    'FontSize',12,...
    'zlim', [0 12000],...
    'ylim', [100 500],...
    'xlim', [0 130]);

% set view
view(ax1,[60 25]);

% Set axes to plot on
axes(ax1);

for aa = 1:5
    % aa = 3
    
    % Assemble X, Y, Z matrices
    [X,Y] = meshgrid(cases2print.rpm{aa},-cases2print.dpdx{aa});
    Z = zeros(length(cases2print.dpdx{aa}),length(cases2print.rpm{aa}));
    for bb=1:length(cases2print.rpm{aa})
        if cases2print.Re_MR{aa,bb}>0
            Z(:,bb) = cases2print.Re_MR{aa,bb}';
        else 
            Z(:,bb) = nan;
        end
    end
    
    % Plot ROP
    if aa<4
        surf(X,Y,Z,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5);
        if aa==1
            text(max(max(X)),max(max(Y)),min(max(Z)),'e-','Color','k','FontSize',12);
        elseif aa==2
            text(max(max(X)),max(max(Y)),min(max(Z)),'e0','Color','k','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'e+','Color','k','FontSize',12);
        end
    else
        surf(X,Y,Z,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r');
        if aa==4
            text(max(max(X)),max(max(Y)),min(max(Z)),'AW','Color','r','FontSize',12);
        else
            text(max(max(X)),max(max(Y)),min(max(Z)),'SW','Color','r','FontSize',12);
        end
    end
end

% set(gca,'zscale', 'log');

zlabel('Re_M_R [-]');

TightFigure(1);

% Print to files
fig_name = 'MP_Re-rpm-dpdx_PAC1';
set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
print(gcf,[fig_path fig_name],'-dsvg'); 






% %% Plot vs. dpdx for 130 rpm and PAC1
% 
% fig = CreateFigure('U_f, U_s, CTR vs. dpdx for 130 rpm and PAC1','Pressure gradient -\Deltap/\Deltax [Pa/m]','Time averaged U_f, U_s  [m/s]', {'lin' 'lin'}, 'DINA5');
% 
% % Select case2print, last gets indices of cell array elements w/ lengths > 1
% cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'PAC1'),:);
% 
% % Select rpm
% bb = 5; % 130 rpm
% 
% % Specify axis colors
% left_color = [.5 .5 0];
% right_color = [0 .5 .5];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% % get(gcf,'defaultAxesColorOrder')
% % % left axis blue [0 0.4470 0.7410]
% % right axis red [0.8500 0.3250 0.0980]
% 
% % Set axis limits
% set(gca,...
%     'xlim', [100 550],...
%     'ylim', [1e-3 2e0]);
% 
% yyaxis left
% 
% for aa = 1:5
%     % aa = 2
% 
%     % Determine sign of plot eccentricity
%     if (cases2print.e(aa)<0 && aa<4)
% %         sign = -1;
%         markersymbol = 'v';
%     elseif (cases2print.e(aa)==0 && aa<4)
% %         sign = 0;
%         markersymbol = 'o';
%     elseif (cases2print.e(aa)>0 && aa<4)
% %         sign = 1;
%         markersymbol = '^';
%     else
%         if aa==4
%              markersymbol = '<';
%         else
%              markersymbol = '>';
%         end
%     end    
%     
%     % Plot fluid superficial velocity vs dpdx
%     plot(-cases2print.dpdx{aa},cases2print.U_f_x{aa,bb},...
%         'LineStyle','-',...
%         'Color',[0 0.4470 0.7410],...
%         'Marker', markersymbol,...
%         'MarkerSize',5,...
%         'MarkerEdgeColor',[0 0.4470 0.7410],...
%         'MarkerFaceColor',[0 0.4470 0.7410]);
%     % 
%     % 'LineWidth',turbmodels{idx_turb,3},...
%     
%     % Determine scale of markers relative to axis
% %     ax = gca;
% %     old_units = get(ax, 'Units');
% %     set(ax, 'Units', 'points');
% %     pos_points = get(ax, 'Position');
% %     set(ax, 'Units', old_units);
% %     scale = max(ylim)/min(pos_points(3:4));
% %     delta = 1*scale;
%     
% 
%     
%     % !!! Annulus marker work, but only for lin plot. for log, somehow use log10 to convert delta 
%     % Use sign, first marker 10, second marker 5
%     
%     % Plot drill pipe as white point
% %     nolegend = plot(-cases2print.dpdx{aa},cases2print.U_f_x{aa,bb}+sign*(delta),...
% %         'LineStyle','none',...
% %         'Marker', 'o',...
% %         'MarkerSize',5,...
% %         'MarkerEdgeColor',[0 0.4470 0.7410],...
% %         'MarkerFaceColor',[0 0.4470 0.7410]);
% %     set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     
%     
%     
%     % Plot solid superficial velocity vs dpdx
%     plot(-cases2print.dpdx{aa},cases2print.U_s_x{aa,bb},...
%         'LineStyle','-',...
%         'Color',[0.8500 0.3250 0.0980],...
%         'Marker', markersymbol,...
%         'MarkerSize',5,...
%         'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
%         'MarkerFaceColor',[0.8500 0.3250 0.0980]);
%     % 
%     % 'LineWidth',turbmodels{idx_turb,3},...
%     
%     % Plot drill pipe as white point
% %     nolegend = plot(-cases2print.dpdx{aa},cases2print.U_s_x{aa,bb}+sign*delta,...
% %         'LineStyle','none',...
% %         'Marker', 'o',...
% %         'MarkerSize',5,...
% %         'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
% %         'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% %     set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     
%     
% end
% 
% set(gca,...
%     'yscale', 'log',...
%     'YColor','k');
% ylabel('\fontsize{11}Time averaged {\color[rgb]{0 0.4470 0.7410}U_f}, {\color[rgb]{0.8500 0.3250 0.0980}U_s} [m/s]','interpreter','tex');
% 
% 
% yyaxis right
% 
% for aa = 1:5
%     % aa = 3
%     
%     % Determine sign of plot eccentricity
%     plotcolor = 'k';
%     if (cases2print.e(aa)<0 && aa<4)
% %         sign = -1;
%         markersymbol = 'v';
%     elseif (cases2print.e(aa)==0 && aa<4)
% %         sign = 0;
%         markersymbol = 'o';
%     elseif (cases2print.e(aa)>0 && aa<4)
% %         sign = 1;
%         markersymbol = '^';
%     else
%         if aa==4
%              markersymbol = '<';
%         else
%              markersymbol = '>';
%         end
%         plotcolor = 'k';
%     end      
%     
%     % Plot fluid superficial velocity vs dpdx
%     plot(-cases2print.dpdx{aa},cases2print.CTR_Q{aa,bb},...
%         'LineStyle','-',...
%         'Color',plotcolor,...
%         'Marker', markersymbol,...
%         'MarkerSize',5,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','k');
%     % 
%     % 'LineWidth',2,...
%     
%     % Determine scale of markers relative to axis
% %     ax = gca;
% %     old_units = get(ax, 'Units');
% %     set(ax, 'Units', 'points');
% %     pos_points = get(ax, 'Position');
% %     set(ax, 'Units', old_units);
% %     scale = max(ylim)/min(pos_points(3:4));
% %     delta = 1*scale;
%                    
%     % Plot drill pipe as white point
% %     nolegend = plot(-cases2print.dpdx{aa},cases2print.CTR_Q{aa,bb}+sign*delta,...
% %         'LineStyle','none',...
% %         'Marker', 'o',...
% %         'MarkerSize',5,...
% %         'MarkerEdgeColor','k',...
% %         'MarkerFaceColor','k');
% %     set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     
% end
% 
% ylabel('CTR [-]');
% 
% set(gca,'yscale', 'log',...
%     'ylim', [1e-3 2e0],...
%     'YColor','k');
% 
% TightFigure(2);
% 
% 
% % Print to files
% fig_name = 'MP_CTR-dpdx_PAC1';
% set(gcf,'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
% print(gcf,[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% % print(gcf,[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
% % print(gcf,[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'
% print(gcf,[fig_path fig_name],'-dsvg');  % '-dpng' '-depsc' 'jpeg'



%% Plot vs time

% % Select case2print, last gets indices of cell array elements w/ lengths > 1
% cases2print = cases(cases.L>0.64 & strcmp(cases.Fluid,'H2O'),:);
% 
% % Select rpm
% bb = 5; % 130 rpm
% 
% 
% CreateFigure('U_f, U_s, CTR vs. time','','', {'lin' 'lin'}, 'DINA5');
% 
% 
% % Whirling simulations have different total flow times than non-whirling
% % due to remeshing intervalls. For time comparisons, the individual flow
% % time intervals of the non-whirling simulation have to be shortened such
% % that they can be plotted together with the local whirling simulations.
% 
% for aa = 1:5
%     % aa = 2
% 
%     % Determine sign of plot eccentricity
%     if (cases2print.e(aa)<0 && aa<4)
% %         sign = -1;
%         markersymbol = 'v';
%     elseif (cases2print.e(aa)==0 && aa<4)
% %         sign = 0;
%         markersymbol = 'o';
%     elseif (cases2print.e(aa)>0 && aa<4)
% %         sign = 1;
%         markersymbol = '^';
%         plotcolor = 'k';
%     else
%         if aa==4
%              markersymbol = '<';
%         else
%              markersymbol = '>';
%              plotcolor = 'r';
%         end
%     end
%     
%     % Plot fluid superficial velocity vs dpdx
%         %yyaxis left
%         subplot(2,1,1); hold on;
%         plot(cases2print.timeseries{aa},cases2print.U_f_x_series{aa,bb}(:,1),...
%             'LineStyle','-',...
%             'Color',plotcolor,...
%             'Marker', '.',...
%             'MarkerSize',10,...
%             'MarkerEdgeColor',plotcolor,...
%             'MarkerFaceColor',plotcolor);
% 
%         % Plot solid superficial velocity vs dpdx
%         subplot(2,1,2); hold on;
%         plot(cases2print.timeseries{aa},cases2print.U_s_x_series{aa,bb}(:,1),...
%             'LineStyle','-',...
%             'Color',plotcolor,...
%             'Marker', '.',...
%             'MarkerSize',10,...
%             'MarkerEdgeColor',plotcolor,...
%             'MarkerFaceColor',plotcolor);
% 
% end
% 
% 
% for ii = 1:2
% subplot(2,1,ii);
% xlabel('Pressure gradient -\Deltap/\Deltax [Pa/m]');
% if ii==1
%     ylabel('(Time averaged) U_f [m/s]');
%     set(gca,'Ylim',[0 1.5]);
% else
%     ylabel('(Time averaged) U_s [m/s]');
%     set(gca,'Ylim',[0 0.05]);
% end
% box on;
% grid on;
% TightFigure(1);
% end
% 
% 
% 
% % set(gca,...
% %     'xlim', [0 round(time(end),0)],...
% %     'ylim', [0 round(U_f(end),2)]);
% 
% 
% yyaxis right
% 
% dpdx = (100:100:500)';
% 
% dt = 1e-3;
% tmax = GetRealFlowTime( 30 );
% idx_low = 20*(1/dt);
% idx_up = tmax*(1/dt);
% time2=dt*[0; repelem((1:1:length(dpdx)-1)',2,1); length(dpdx)]*idx_up;
% dpdx2 = repelem(dpdx,2,1);
% 
% for ii = 1:2
%     subplot(2,1,ii);
%     plot(time2,dpdx2);
%     set(gca,...
%         'xlim', [0 round(time(end),0)],...
%         'ylim', [0 dpdx2(end)]);
%     ylabel('-\Deltap/\Deltax [Pa/m]');
% end
% 
% TightFigure(2);
