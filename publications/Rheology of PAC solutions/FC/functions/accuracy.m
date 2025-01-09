
% set(gca,'XLim',[1e-3 1e4])
% set(gca,'XLim',[1e-2 1e3])
SR = logspace(-3,4);

% Rheometer
r_i = 13.33/1000; 
r_o = 14.46/1000;
H = 39.99/1000;

% Fluid sample
rho = 1000;


%% Lower limit based on minimum torque
% based on Irgens (2014)

% Accuracy +/- 5 %, based on pptx in Email Jan Schaffer (18.09.2017)
M  = 5e-9;
factor = M./(2.*pi.*r_i.^2.*(H+r_i./3));
eta = factor./SR;
limit_low_5 = plot(SR,eta,'k-');
% text(2e-3,0.02,'±5%','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Color','k'); %'interpreter','latex'
% annotation('textarrow',[0.15 0.1],[0.2 0.2],'String','','Color','k');
text(4e0,0.02,'±5%','HorizontalAlignment','right','FontSize',12,'FontWeight','bold','Color','k'); %'interpreter','latex'
annotation('textarrow',[0.61 0.55],[0.23 0.23],'String','','Color','k');


% Accuracy +/- 1 %, based on Email Jan Schaffer (18.09.2017)
M  = 10e-6;
factor = M./(2.*pi.*r_i.^2.*(H+r_i./3));
eta = factor./SR;
limit_low_1 = plot(SR,eta,'k--');
text(2.5e1,0.02,'±1%','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Color','k'); %'interpreter','latex'
annotation('textarrow',[0.635 0.69],[0.23 0.23],'String','','Color','k');


% Different formula from Chhabra (2008), probably no bottom torque
% eta = M./(4.*pi.*H.*SR.*r_i.*r_o);
% plot(SR,eta,'Color','k');

annotation('textarrow',[0.2 0.09],[0.18 0.12],'String','>±5%','Color','k','FontSize',12);

%% Upper limit based on Taylor vortices / turbulence

% Rodd (2005)
eta = rho*SR.^2.*sqrt((r_o-r_i)^7)./(1700.*r_i);
% limit_high_taylor = plot(SR,eta,'k');

% Barnes (2000)
ta_cr = 1700*(1+r_o/10/r_i);
factor = sqrt(7.036521384*10.^(-11));
eta = factor.*SR;
limit_high_taylor = plot(SR,eta,'Color','k');

annotation('textarrow',[0.885 0.986],[0.18 0.12],'String','Taylor vortices','Color','k','FontSize',12);



%% Legend

% Toggle legend display status
set(get(get(limit_low_5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(limit_low_1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(limit_high_taylor,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');