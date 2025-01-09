%% Get data

% run RUN2.m with
% resultspath = 'E:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - 2D\archive\170914_2D_C_expVelProf_slip_Cross_SR=1\exports';
for i=1:2 % i =2
    y_nof{i} = CFD_asis.Y{i,1};
    x_nof{i} = CFD_asis.X{i,1}-CFD_asis.X{i,1}(1);
    y_nof{i}(any(isnan(y_nof{i}),2),:)=[];
    x_nof{i}(any(isnan(x_nof{i}),2),:)=[];
end
    
% run RUN2.m with
% resultspath = 'O:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - 2D\archive\171019_2D_C_expVelProf_slip_Cross_SR=1_VM_PressGrad_Lift\exports';
% ! Comment out clear all; !
for i=1:2
    y_of{i} = CFD_asis.Y{i,1};
    x_of{i} = CFD_asis.X{i,1}-CFD_asis.X{i,1}(1);
    y_of{i}(any(isnan(y_of{i}),2),:)=[];
    x_of{i}(any(isnan(x_of{i}),2),:)=[];
end

%% Plot

close all;
fig_OtherForces = CreateFigure( 'Relevance of other forces', 'Particle y-coordinate (= channel y-coordinate) [m]', 'Particle x-coordinate ratio x_o_f/x [-]' );
set(gcf,'Position',[1 1 21 14.8/2]);
hold on;
plotcolor = colorlist{1};

y = linspace(0.02,0)';

for i=1:2
    plot(y,interp1(y_of{i},x_of{i},y,'pchip') ./ interp1(y_nof{i},x_nof{i},y,'pchip'),'Color',plotcolor,'LineStyle','-');
end

set(gca,...
    'FontSize',sub_fontsize,...
    'XLim',[0 0.02],...
    'XDir','reverse');

text(0.016,1,{'\color{blue}d_p = 0.00116 m'},'HorizontalAlignment' ,'left','VerticalAlignment','Bottom','FontSize',14);
text(0.012,0.6,{'\color{blue}d_p = 0.002 m'},'HorizontalAlignment' ,'left','VerticalAlignment','Top','FontSize',14);

% Expand Axes to Fill Figure (--> Minimum white space)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
factor_hor = 0; %.003;
factor_ver = 0; %.01;
left = outerpos(1) + ti(1) + factor_hor;
bottom = outerpos(2) + ti(2) + factor_ver;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*factor_hor;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*factor_ver;
ax.Position = [left bottom ax_width ax_height];


%% Save figure
figurepath = '\\home.ansatt.ntnu.no\alexabus\Documents\AdWell\8 - Publications & Conferences\2017-xx - An Eulerian-Lagrangian CFD study of a particle settling in an orthogonal GNF shear flow\Figure';

savefig(fig_OtherForces, [figurepath '\OtherForces.fig']);
print(fig_OtherForces,[figurepath '\OtherForces'],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_OtherForces,[figurepath '\OtherForces'],'-dmeta');  % '-dpng' '-depsc' 'jpeg'

set(gcf, 'color', 'none');
set(findall(gcf,'type','axes'), 'color', 'none');
export_fig OtherForces.png;
movefile('OtherForces.png', figurepath);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');

    