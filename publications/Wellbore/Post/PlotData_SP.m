close all;
clear all;
clc;

% Set all relevant paths
addpath(genpath('C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic'));
addpath([pwd '\utilities']);
fig_path = 'M:\Documents\AdWell\8 - Publications & Conferences\2018-xx - Wellbore\Figures\';
filepath_simdata = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive';
filepath_metadata = [filepath_simdata  '\' 'metadata2.xlsx'];

% Read data
 load('wellboredata.mat');

 % Create subset of data to be plotted
% cases2print = cases(strcmp(cases.Job,'Single') ...
%      & contains(cases.Purpose,'pressure sweep')==1 ...
%      ,:);
   
  
cases2print = cases(cases.d_s==0 ...
     & contains(cases.Purpose,'pressure sweep')==1 ...
     ,:); 

 cases2print = sortrows(cases2print,{'e','rpm','Fluid','r','cells'},{'descend','ascend','ascend','descend','ascend'});
 
% Define colorband
N = 5; 
ColorSet = ColorBand(N);
    
 %% Plot 4000 vs 32000 cells
 
 % Create plot
 fig_handle(1) = CreateFigure('\Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_s_f [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
 
 
 legendentries = cell(height(cases2print),0);
 for aa = 1:height(cases2print)
     if (strcmp(cases2print.Fluid{aa},'PAC1') || cases2print.r(aa)==0)
         
     else    
         if cases2print.cells(aa)==4000
             linetype = '--';
         else
             linetype = '-';
             set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
         end
         legendentries{aa,1} = [ 'e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];
         plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'Marker','o');
     end
 end
 
 
% plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'Marker','o');
 
 
% Adjust axis
set(gca,...
    'FontSize',12,...
    'xlim', [0 2.5],...
    'ylim', [0 3500]);

% Legend 
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','northwest')



%% Plot CFD vs UiS exp data - H2O

% Create plot
fig_handle(2) = CreateFigure('\Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_f_,_x [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
legendentries = cell(1,0);
  
% Concentric case
aa=2;

% CFD
linetype = '-';
plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{1,1} = [ 'CFD, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Friction factors
vel = linspace(0.1,2)';
Re = cases2print.rho_f(aa).*vel.*cases2print.d_h(aa)./cases2print.K_PL_f(aa);

% Blasius
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=0.316.*Re.^-0.25/4;
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
plot(vel,dpdx,'LineStyle',linetype);
legendentries{2,1} = 'Blasius';
% 
% % Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% f_Morrison = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
% dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{3,1} = [ 'Morrison, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Haaland
linetype = '-.'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(aa)).^1.11+(6.9./Re))).^2);
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
plot(vel,dpdx,'LineStyle',linetype);
legendentries{3,1} = 'Haaland';


% Eccentric case
aa=6;

% CFD
linetype = '-';
plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{4,1} = [ 'CFD, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Blasius
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=0.316.*Re.^-0.25/4;
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_f(aa)  );
plot(vel,dpdx,'LineStyle',linetype);
legendentries{5,1} = 'Blasius w/ Haciislamoglu & Cartalos (1994) correction';

% % Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% f_Morrison = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
% dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
% dpdx = EccentricPressureCorrection( dpdx, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_f(aa)  );
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{7,1} = [ 'Morrison, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Haaland
linetype = '-.'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(aa)).^1.11+(6.9./Re))).^2);
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_f(aa)  );
plot(vel,dpdx,'LineStyle',linetype);
legendentries{6,1} = 'Haaland w/ Haciislamoglu & Cartalos (1994) correction';

% Experimental data of UiS
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','A9:B21');
scatter(expdata(:,1),expdata(:,2));
% scatter(cases2print.U_f_x_exp{aa,:},-cases2print.dpdx_exp{aa,:});
legendentries{7,1} = [ 'Exp. UiS, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ];

% Experimental data of UiS, corrected with L = 1.52
scatter(expdata(:,1),expdata(:,2)/1.52,'k');
legendentries{8,1} = 'Exp. UiS, corrected w/ L = 1.52 m' ;

% Experimental data of UiS, corrected with L = 1.52
scatter(expdata(:,1),expdata(:,2)/1.52,'k');
legendentries{8,1} = 'Exp. UiS, corrected w/ L = 1.52 m' ;

% Adjust axis
set(gca,...
        'FontSize',12,...
    'xlim', [0 1.5],...
    'ylim', [0 3000]);

% Legend
TightFigure(1);
legend(legendentries,'Location','northwest');

% 
% % Eccentric case CFD case 
% ID = '180926_5';
% 
% 
% % CFD
% % sweep of fluid superficial velocities with transient simulations
% idx = find(strcmp('180926_5',cases.ID(:)));
% linetype = '-';
% plot(cases.U_f_x{idx,:},-cases.dpdx{idx,:},'LineStyle',linetype,'Marker','o');
% legendentries{4,1} = [ 'CFD, e = ' num2str(cases.e(aa)) ', rpm = ' num2str(cases.rpm(aa)) ', ' cases.Fluid{aa} ', r =' num2str(cases.r(aa)) ', cells = ' num2str(cases.cells(aa)) ];



%% Plot CFD vs UiS exp data - PAC1

% Create plot
fig_handle(3) = CreateFigure('PAC1 - \Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_s_f [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
legendentries = cell(1,0);
  

% Concentric case
aa=3;

% CFD
linetype = '-';
plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{1,1} = [ 'CFD, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% Friction factors
vel = linspace(0.1,2)';
Re = cases2print.rho_f(aa).*vel.^(2-cases2print.n_PL_f(aa)).*cases2print.d_h(aa).^cases2print.n_PL_f(aa)./(12.^(cases2print.n_PL_f(aa)-1).*cases2print.K_PL_f(aa).*((2.*cases2print.n_PL_f(aa)+1)./(3.*cases2print.n_PL_f(aa))).^cases2print.n_PL_f(aa));

% % Critical Re as f(n_PL_f)
% figure; hold on;
% % Dodge & Metzner (1959)
% plot(n,3250-1150*n);
% % Ryan and Johnson (1959)
% plot(n,6464*n./(3*n+1).^2.*(2+n).^((2+n)./(1+n)));
% % Mishra and Tripathi (1975)
% plot(n,2100*(4*n+2).*(5*n+3)./(3*(3*n+1).^2));
% legend('Dodge & Metzner (1959)',...
%     'Ryan and Johnson (1959)',...
%     'Mishra and Tripathi (1975)');
Re_cr = 2100*(4*cases2print.n_PL_f(aa)+2).*(5*cases2print.n_PL_f(aa)+3)./(3*(3*cases2print.n_PL_f(aa)+1).^2);
vel_cr = (Re_cr.*(12.^(cases2print.n_PL_f(aa)-1).*cases2print.K_PL_f(aa).*((2.*cases2print.n_PL_f(aa)+1)./(3.*cases2print.n_PL_f(aa))).^cases2print.n_PL_f(aa))./(cases2print.rho_f(aa).*cases2print.d_h(aa).^cases2print.n_PL_f(aa))).^(1./(2-cases2print.n_PL_f(aa)));

% PL, turbulent regime (Irvine 1988)
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=zeros(length(Re),1);
for ii = 1:length(Re)
   if Re(ii)<Re_cr
       f(ii)=24/Re(ii);
   else
       f(ii)=((((2.^(cases2print.n_PL_f(aa)+4))./(7.^(7.*cases2print.n_PL_f(aa)))).*(4.*cases2print.n_PL_f(aa)./(3.*cases2print.n_PL_f(aa)+1)).^(3.*cases2print.n_PL_f(aa).^2))./Re(ii)).^(1./(3.*cases2print.n_PL_f(aa)+1));
   end
end
f_cr = 24/Re_cr;
f_cr = ((((2.^(cases2print.n_PL_f(aa)+4))./(7.^(7.*cases2print.n_PL_f(aa)))).*(4.*cases2print.n_PL_f(aa)./(3.*cases2print.n_PL_f(aa)+1)).^(3.*cases2print.n_PL_f(aa).^2))./Re_cr).^(1./(3.*cases2print.n_PL_f(aa)+1));
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
plot(vel,dpdx,'LineStyle',linetype);
legendentries{2,1} = '(Irvine 1988)';

% figure ,hold on;
% plot(Re,f)
% set(gca,'Xscale','log','Yscale','log');

% PL, turbulent regime (Dodge & Metzner 1959)
% linetype = '-.'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% for ii = 1:length(Re)
%    if f(ii)>f_cr
%        Re(ii)=24/f(ii);
%    else
%        Re(ii) = (10.^((1./sqrt(f(ii))+0.4./(cases2print.n_PL_f(aa).^1.2)).*(cases2print.n_PL_f(aa).^0.75)./4))./(f(ii).^((2-cases2print.n_PL_f(aa))./2));
%    end
% end
% vel = (Re.*(12.^(cases2print.n_PL_f(aa)-1).*cases2print.K_PL_f(aa).*((2.*cases2print.n_PL_f(aa)+1)./(3.*cases2print.n_PL_f(aa))).^cases2print.n_PL_f(aa))./(cases2print.rho_f(aa).*cases2print.d_h(aa).^cases2print.n_PL_f(aa))).^(1./(2-cases2print.n_PL_f(aa)));
% dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{6,1} = '(Dodge & Metzner 1959) w/ Haciislamoglu & Cartalos (1994) correction';


% Eccentric case
aa=9;

% CFD
linetype = '-';
plot(cases2print.U_f_x{aa,1},-cases2print.dpdx{aa,:},'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{4,1} = [ 'CFD, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ', cells = ' num2str(cases2print.cells(aa)) ];

% PL, turbulent regime (Irvine 1988)
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=zeros(length(Re),1);
for ii = 1:length(Re)
   if Re(ii)<Re_cr
       f(ii)=24/Re(ii);
   else
       f(ii)=((((2.^(cases2print.n_PL_f(aa)+4))./(7.^(7.*cases2print.n_PL_f(aa)))).*(4.*cases2print.n_PL_f(aa)./(3.*cases2print.n_PL_f(aa)+1)).^(3.*cases2print.n_PL_f(aa).^2))./Re(ii)).^(1./(3.*cases2print.n_PL_f(aa)+1));
   end
end
f_cr = 24/Re_cr;
f_cr = ((((2.^(cases2print.n_PL_f(aa)+4))./(7.^(7.*cases2print.n_PL_f(aa)))).*(4.*cases2print.n_PL_f(aa)./(3.*cases2print.n_PL_f(aa)+1)).^(3.*cases2print.n_PL_f(aa).^2))./Re_cr).^(1./(3.*cases2print.n_PL_f(aa)+1));
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_f(aa)  );
plot(vel,dpdx,'LineStyle',linetype);
legendentries{5,1} = '(Irvine 1988) w/ Haciislamoglu & Cartalos (1994) correction';

% PL, turbulent regime (Dodge & Metzner 1959)
% linetype = '-.'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% for ii = 1:length(Re)
%    if f(ii)>f_cr
%        Re(ii)=24/f(ii);
%    else
%        Re(ii) = (10.^((1./sqrt(f(ii))+0.4./(cases2print.n_PL_f(aa).^1.2)).*(cases2print.n_PL_f(aa).^0.75)./4))./(f(ii).^((2-cases2print.n_PL_f(aa))./2));
%    end
% end
% vel = (Re.*(12.^(cases2print.n_PL_f(aa)-1).*cases2print.K_PL_f(aa).*((2.*cases2print.n_PL_f(aa)+1)./(3.*cases2print.n_PL_f(aa))).^cases2print.n_PL_f(aa))./(cases2print.rho_f(aa).*cases2print.d_h(aa).^cases2print.n_PL_f(aa))).^(1./(2-cases2print.n_PL_f(aa)));
% dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
% dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, cases2print.d_i(aa), cases2print.d_o(aa), -cases2print.e(aa), cases2print.n_PL_f(aa)  );
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{6,1} = '(Dodge & Metzner 1959) w/ Haciislamoglu & Cartalos (1994) correction';

% Experimental data of UiS
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','D9:E30');
scatter(expdata(:,1),expdata(:,2));
% scatter(cases2print.U_f_x_exp{aa,:},-cases2print.dpdx_exp{aa,:});
legendentries{7,1} = [ 'Exp. UiS, e = ' num2str(cases2print.e(aa)) ', rpm = ' num2str(cases2print.rpm{aa}) ', ' cases2print.Fluid{aa} ', r =' num2str(cases2print.r(aa)) ];

% Experimental data of UiS, corrected with L = 1.52
scatter(expdata(:,1),expdata(:,2)/1.52,'k');
legendentries{8,1} = 'Exp. UiS, corrected w/ L = 1.52 m' ;

% Adjust axis
set(gca,...
        'FontSize',12,...
    'xlim', [0 2.2],...
    'ylim', [0 4500]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','northwest');



%% Plot CFD vs UiS exp data - H2O

% Create plot
fig_handle(4) = CreateFigure('\Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_f_,_x [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
fig_handle(5) = CreateFigure('f vs. Re_MR','Metzner-Reed Reynolds number Re_M_R [-]', 'Friction factor f [-]', {'log' 'log'}, 'DINA5');
legendentries = cell(1,0);

% Set colorband
figure(fig_handle(4));
set(gca,'ColorOrder', ColorSet);
figure(fig_handle(5));
set(gca,'ColorOrder', ColorSet);

% Concentric case

% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.e==0 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'dpdx')...
    & cellfun(@(x) length(x) > 1,cases.dpdx) ,:);


% CFD
linetype = '-';
figure(fig_handle(4));
plot(case2print.U_f_x{1,1},-case2print.dpdx{:},'LineStyle',linetype,'LineWidth',2,'Marker','o');
figure(fig_handle(5));
plot(case2print.Re_MR{1,1},-case2print.f{1,1},'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{1,1} = [ 'CFD, e = ' num2str(case2print.e) ', r =' num2str(cases2print.r(aa)) ' m' ];


% Friction factors
vel = linspace(0.1,2.5)';
Re = case2print.rho_f.*vel.*case2print.d_h./case2print.K_PL_f;

% Blasius
linetype = '--';
f=0.316.*Re.^-0.25/4;
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(vel,dpdx,'LineStyle',linetype);
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(Re,f,'LineStyle',linetype);
legendentries{2,1} =[ 'Blasius (1912), e = ' num2str(case2print.e) ', r = 0 m' ];

% plot(Re,24./Re,'k--')
% Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% f_Morrison = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
% dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{3,1} = [ 'Morrison (2013), e = ' num2str(case2print.e) ', r = 0 m' ];

% Haaland
linetype = '-.';
f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(aa)).^1.11+(6.9./Re))).^2);
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(vel,dpdx,'LineStyle',linetype);
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(Re,f,'LineStyle',linetype);
legendentries{4,1} = [ 'Haaland (1983), e = ' num2str(case2print.e) ', r = ' num2str(cases2print.r(aa)) ' m' ];

% Eccentric case

% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.e==-0.95 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'dpdx')...
    & cellfun(@(x) length(x) > 1,cases.dpdx) ,:);

case2print2 = cases(cases.e==-0.95 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'U_f')...
    & cellfun(@(x) length(x) > 1,cases.U_f_x) ,:);


x = [case2print2.U_f_x{1,1}(case2print2.U_f_x{1,1} < min(case2print.U_f_x{1,1})); case2print.U_f_x{1,1}];
y = -[case2print2.dpdx{1,1}(case2print2.U_f_x{1,1} < min(case2print.U_f_x{1,1})); case2print.dpdx{1,1}];
y(1)=117; y(2)=232; x(3)=[]; y(3)=[]; % ScanIT from Figure 6 of 1st submission


% CFD
linetype = '-';
figure(fig_handle(4));
plot(x,y, 'LineStyle',linetype,'LineWidth',2,'Marker','o');
%plot(case2print.U_f_x{1,1},case2print.dpdx{1,1}, 'LineStyle',linetype,'LineWidth',2,'Marker','o');

legendentries{5,1} = [ 'CFD, e = ' num2str(cases2print.e(aa)) ', r =' num2str(cases2print.r(aa)) ' m' ];
figure(fig_handle(5));
x = [case2print2.Re_MR{1,1}(case2print2.Re_MR{1,1} < min(case2print.Re_MR{1,1})); case2print.Re_MR{1,1}];
y = -[case2print2.f{1,1}(case2print2.Re_MR{1,1} < min(case2print.Re_MR{1,1})); case2print.f{1,1}];
%plot(x,y,'LineStyle',linetype,'LineWidth',2,'Marker','o');
plot(case2print.Re_MR{1,1},case2print.f{1,1},'LineStyle',linetype,'LineWidth',2,'Marker','o');

clear x y case2print2

% Blasius
linetype = '--';
f=0.316.*Re.^-0.25/4;
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(vel,dpdx,'LineStyle',linetype);
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(Re,dpdx./(case2print.rho_f.*vel.^2./(case2print.d_h./2)),'LineStyle',linetype);
legendentries{6,1} = ['Blasius (1912), e = ' num2str(cases2print.e(aa)) ', r = 0 m' ];

% % Morrison
% linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
% f_Morrison = (0.0076*(3170./Re).^0.165)./(1+(3170./Re).^7)+16./Re;
% dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
% dpdx = EccentricPressureCorrection( dpdx, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');
% legendentries{7,1} = [ 'Morrison, e = ' num2str(case2print.e) ', rpm = ' num2str(case2print.rpm) ', ' case2print.Fluid{aa} ', r =' num2str(case2print.r) ', cells = ' num2str(case2print.cells) ];

% Haaland
linetype = '-.';
f = 1./(4.* (-1.8.*log10((0.00003./3.7/case2print.d_h).^1.11+(6.9./Re))).^2);
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(vel,dpdx,'LineStyle',linetype);
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(Re,dpdx./(case2print.rho_f.*vel.^2./(case2print.d_h./2)),'LineStyle',linetype);
legendentries{7,1} = ['Haaland (1983), e = ' num2str(cases2print.e(aa)) ', r = ' num2str(cases2print.r(aa)) ' m' ];

% Experimental data of UiS
expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','A9:B21');
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(expdata(:,1),expdata(:,2));
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(case2print.rho_f.*expdata(:,1).*case2print.d_h./case2print.K_PL_f,expdata(:,2)./(case2print.rho_f.*expdata(:,1).^2./(case2print.d_h./2)));
legendentries{8,1} = [ 'Khatibi et al. (2018a), e = ' num2str(cases2print.e(aa)) ];

expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','C58:D139');
figure(fig_handle(4));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(expdata(:,1),expdata(:,2),'.');
figure(fig_handle(5));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(case2print.rho_f.*expdata(:,1).*case2print.d_h./case2print.K_PL_f,expdata(:,2)./(case2print.rho_f.*expdata(:,1).^2./(case2print.d_h./2)),'.');
legendentries{9,1} = [ 'Khatibi et al. (2018b), e = ' num2str(cases2print.e(aa)) ];


% dpdx-U
figure(fig_handle(4));

% Adjust axis
set(gca,...
        'FontSize',12,...
    'xlim', [0 2.5],...
    'ylim', [0 3500]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','northwest');

% Print to files
fig_name = 'SP_dpdx-U';
set(fig_handle(4),'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(fig_handle(4),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(4),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(4),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'


% f-Re
figure(fig_handle(5));

% Adjust axis
set(gca,...
    'xlim', [4000 40000],...
    'ylim', [0.002 0.02]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','southwest');

% Print to files
fig_name = 'SP_f-Re';
set(fig_handle(5),'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(fig_handle(5),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(5),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(5),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'


%% Plot CFD vs UiS exp data - H2O - Whirling motion

% Create plot
fig_handle(6) = CreateFigure('\Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_f_,_x [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
fig_handle(7) = CreateFigure('f vs. Re_MR','Metzner-Reed Reynolds number Re_M_R [-]', 'Friction factor f [-]', {'log' 'log'}, 'DINA5');
legendentries = cell(1,0);

% Set colorband
figure(fig_handle(6));
set(gca,'ColorOrder', ColorSet);
figure(fig_handle(7));
set(gca,'ColorOrder', ColorSet);

% Concentric case as benchmark

% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.e==0 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'dpdx')...
    & cellfun(@(x) length(x) > 1,cases.dpdx) ,:);

% Friction factors
vel = linspace(0.1,2.5)';
Re = case2print.rho_f.*vel.*case2print.d_h./case2print.K_PL_f;

% Blasius
linetype = '--'; 
f=0.316.*Re.^-0.25/4;
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
figure(fig_handle(6));
nolegend = plot(vel,dpdx,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
nolegend = plot(Re,f,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% Haaland
linetype = '-.';
f = 1./(4.* (-1.8.*log10((0.00003./3.7/cases2print.d_h(aa)).^1.11+(6.9./Re))).^2);
dpdx = f.*cases2print.rho_f(aa).*vel.^2./(cases2print.d_h(aa)./2);
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(vel,dpdx,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(Re,f,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% Eccentric cases with 0 rpm

% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.e==-0.95 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'dpdx')...
    & cellfun(@(x) length(x) > 1,cases.dpdx) ,:);

case2print2 = cases(cases.e==-0.95 ...
    & cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.flow_spec,'U_f')...
    & cellfun(@(x) length(x) > 1,cases.U_f_x) ,:);


x = [case2print2.U_f_x{1,1}(case2print2.U_f_x{1,1} < min(case2print.U_f_x{1,1})); case2print.U_f_x{1,1}];
y = -[case2print2.dpdx{1,1}(case2print2.U_f_x{1,1} < min(case2print.U_f_x{1,1})); case2print.dpdx{1,1}];
y(1)=117; y(2)=232; x(3)=[]; y(3)=[]; % ScanIT from Figure 6 of 1st submission

% CFD
linetype = '-';
figure(fig_handle(6));
nolegend = plot(x,y, 'LineStyle',linetype,'LineWidth',2,'Marker','o');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
x = [case2print2.Re_MR{1,1}(case2print2.Re_MR{1,1} < min(case2print.Re_MR{1,1})); case2print.Re_MR{1,1}];
y = -[case2print2.f{1,1}(case2print2.Re_MR{1,1} < min(case2print.Re_MR{1,1})); case2print.f{1,1}];
nolegend = plot(x,y,'LineStyle',linetype,'LineWidth',2,'Marker','o');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

clear x y case2print2

% Blasius
linetype = '--';
f=0.316.*Re.^-0.25/4;
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(vel,dpdx,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(Re,dpdx./(case2print.rho_f.*vel.^2./(case2print.d_h./2)),'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Haaland
linetype = '-.'; 
f = 1./(4.* (-1.8.*log10((0.00003./3.7/case2print.d_h).^1.11+(6.9./Re))).^2);
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, 2100, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(vel,dpdx,'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = plot(Re,dpdx./(case2print.rho_f.*vel.^2./(case2print.d_h./2)),'LineStyle',linetype);
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Experimental data of UiS - 0 rpm
expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','A9:B21');
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(expdata(:,1),expdata(:,2));
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(case2print.rho_f.*expdata(:,1).*case2print.d_h./case2print.K_PL_f,expdata(:,2)./(case2print.rho_f.*expdata(:,1).^2./(case2print.d_h./2)));
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','C60:D141');
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(expdata(:,1),expdata(:,2),'.');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(case2print.rho_f.*expdata(:,1).*case2print.d_h./case2print.K_PL_f,expdata(:,2)./(case2print.rho_f.*expdata(:,1).^2./(case2print.d_h./2)),'.');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','H44:R51');
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter([expdata(1,1) expdata(5,1)],[expdata(1,end) expdata(5,end)],'d');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(case2print.rho_f.*[expdata(1,1) expdata(5,1)].*case2print.d_h./case2print.K_PL_f,[expdata(1,end) expdata(5,end)]./(case2print.rho_f.*[expdata(1,1) expdata(5,1)].^2./(case2print.d_h./2)),'d');
set(get(get( nolegend,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Eccentric cases with 100 rpm

% CFD
% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.L==0.04 ...
    & cases.e<0 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 100') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.CurvCorr,'Yes') ...
    & strcmp(cases.IntTrans,''),:);
linetype = '-';
figure(fig_handle(6));
plot([case2print.U_f_x{:}],-[case2print.dpdx{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
figure(fig_handle(7));
plot([case2print.Re_MR{:}],-[case2print.f{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{1,1} = [ 'CFD, ' num2str(case2print.rpm{1}) ' rpm' ];

% Experimental data of UiS
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter([expdata(2,1) expdata(6,1)],[expdata(2,end) expdata(6,end)],'d');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(case2print.rho_f.*[expdata(2,1) expdata(6,1)]'.*case2print.d_h./case2print.K_PL_f,[expdata(2,end) expdata(6,end)]'./(case2print.rho_f.*[expdata(2,1) expdata(6,1)]'.^2./(case2print.d_h./2)),'d');
legendentries{2,1} = [ 'Khatibi et al. (2018b) Figure 4 & 5, ' num2str(case2print.rpm{1}) ' rpm' ];



% Eccentric cases with 200 rpm

% CFD
% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 200') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.CurvCorr,'Yes'),:);
linetype = '-';
figure(fig_handle(6));
plot([case2print.U_f_x{:}],-[case2print.dpdx{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
figure(fig_handle(7));
plot([case2print.Re_MR{:}],-[case2print.f{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{3,1} = [ 'CFD, ' num2str(case2print.rpm{1}) ' rpm' ];

% Experimental data of UiS
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter([expdata(3,1) expdata(7,1)],[expdata(3,end) expdata(7,end)],'d');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(case2print.rho_f.*[expdata(3,1) expdata(7,1)]'.*case2print.d_h./case2print.K_PL_f,[expdata(3,end) expdata(7,end)]'./(case2print.rho_f.*[expdata(3,1) expdata(7,1)]'.^2./(case2print.d_h./2)),'d');
legendentries{4,1} = [ 'Khatibi et al. (2018b) Figure 4 & 5, ' num2str(case2print.rpm{1}) ' rpm' ];



% Eccentric cases with 300 rpm

% CFD
% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.L==0.04 ...
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 300') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'H2O') ...
    & strcmp(cases.CurvCorr,'Yes'),:);
linetype = '-';
figure(fig_handle(6));
plot([case2print.U_f_x{:}],-[case2print.dpdx{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
figure(fig_handle(7));
plot([case2print.Re_MR{:}],-[case2print.f{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');
legendentries{5,1} = [ 'CFD, ' num2str(case2print.rpm{1}) ' rpm' ];

% Experimental data of UiS - 300 rpm
figure(fig_handle(6));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter([expdata(4,1) expdata(8,1)],[expdata(4,end) expdata(8,end)],'d');
figure(fig_handle(7));
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
nolegend = scatter(case2print.rho_f.*[expdata(4,1) expdata(8,1)]'.*case2print.d_h./case2print.K_PL_f,[expdata(4,end) expdata(8,end)]'./(case2print.rho_f.*[expdata(4,1) expdata(8,1)]'.^2./(case2print.d_h./2)),'d');
legendentries{6,1} = [ 'Khatibi et al. (2018b) Figure 4 & 5, ' num2str(case2print.rpm{1}) ' rpm' ];



% dpdx-U
figure(fig_handle(6));

% Adjust axis
set(gca,...
        'FontSize',12,...
    'xlim', [0.32 0.52],...
    'ylim', [100 500]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','northwest');

% Print to files
fig_name = 'SP_dpdx-U_zoom_whirl';
set(fig_handle(6),'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(fig_handle(6),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(6),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(6),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'


% f-Re
figure(fig_handle(7));

% Adjust axis
set(gca,...
    'xlim', [4000 8000],...
    'ylim', [0.005 0.02]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','southwest');

% Print to files
fig_name = 'SP_f-Re_zoom_whirl';
set(fig_handle(7),'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(fig_handle(7),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(7),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(7),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'


%% Plot CFD vs UiS exp data - PAC1 

% Create plot
fig_handle(8) = CreateFigure('PAC1 - \Deltap/\Deltax vs. U_s_f','Superficial fluid velocity U_s_f [m/s]', 'Pressure gradient -\Deltap/\Deltax [Pa/m]', {'lin' 'lin'}, 'DINA5');
legendentries = cell(1,0);
 
% Set colorband
figure(fig_handle(8));
set(gca,'ColorOrder', ColorSet);
% figure(fig_handle(9));
% set(gca,'ColorOrder', ColorSet);


% CFD, concentric case
% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.L==0.04 ...
    & cases.e==-0 ...    
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'PAC1') ...
    & strcmp(cases.CurvCorr,'Yes') ...
    & strcmp(cases.IntTrans,''),:);
linetype = '-';
%plot([case2print.U_f_x{1,1}],-[case2print.dpdx{1,1}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');

% ScanIT Fig7 of 1st submission
x =    [0.0246	0.344	0.529	0.689	0.889	1.1	1.27	1.44]
y = [545	844	1.23E+03	1.59E+03	2.17E+03	2.88E+03	3.61E+03	4.43E+03]
plot(x,y, 'LineStyle',linetype,'LineWidth',2,'Marker','o');

legendentries{1,1} = [ 'CFD, e = ' num2str(case2print.e) ', r = ' num2str(case2print.r) ' m' ];


% Friction factors
vel = linspace(0.1,2.5)';
Re = case2print.rho_f.*vel.^(2-case2print.n_PL_f).*case2print.d_h.^case2print.n_PL_f./(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f);

% % Critical Re as f(n_PL_f)
% figure; hold on;
% % Dodge & Metzner (1959)
% plot(n,3250-1150*n);
% % Ryan and Johnson (1959)
% plot(n,6464*n./(3*n+1).^2.*(2+n).^((2+n)./(1+n)));
% % Mishra and Tripathi (1975)
% plot(n,2100*(4*n+2).*(5*n+3)./(3*(3*n+1).^2));
% legend('Dodge & Metzner (1959)',...
%     'Ryan and Johnson (1959)',...
%     'Mishra and Tripathi (1975)');
Re_cr = 2100*(4*case2print.n_PL_f+2).*(5*case2print.n_PL_f+3)./(3*(3*case2print.n_PL_f+1).^2);
vel_cr = (Re_cr.*(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f)./(case2print.rho_f.*case2print.d_h.^case2print.n_PL_f)).^(1./(2-case2print.n_PL_f));


% PL, turbulent regime (Irvine 1988)
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=zeros(length(Re),1);
for ii = 1:length(Re)
   if Re(ii)<Re_cr
       f(ii)=24/Re(ii);
   else
       f(ii)=((((2.^(case2print.n_PL_f+4))./(7.^(7.*case2print.n_PL_f))).*(4.*case2print.n_PL_f./(3.*case2print.n_PL_f+1)).^(3.*case2print.n_PL_f.^2))./Re(ii)).^(1./(3.*case2print.n_PL_f+1));
   end
end
f_cr = 24/Re_cr;
f_cr = ((((2.^(case2print.n_PL_f+4))./(7.^(7.*case2print.n_PL_f))).*(4.*case2print.n_PL_f./(3.*case2print.n_PL_f+1)).^(3.*case2print.n_PL_f.^2))./Re_cr).^(1./(3.*case2print.n_PL_f+1));
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
plot(vel,dpdx,'LineStyle',linetype);
legendentries{2,1} = ['Irvine (1988), e = ' num2str(case2print.e) ', r = 0 m' ];


% PL, turbulent regime (Dodge & Metzner 1959), implicit in f hence specify
% existing f >> Re >> vel
linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
for ii = 1:length(Re)
   if f(ii)>f_cr
       Re(ii)=24/f(ii);
   else
       Re(ii) = (10.^((1./sqrt(f(ii))+0.4./(case2print.n_PL_f.^1.2)).*(case2print.n_PL_f.^0.75)./4))./(f(ii).^((2-case2print.n_PL_f)./2));
   end
end
vel = (Re.*(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f)./(case2print.rho_f.*case2print.d_h.^case2print.n_PL_f)).^(1./(2-case2print.n_PL_f));
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
plot(vel,dpdx,'LineStyle',linetype);
legendentries{3,1} = ['Dodge & Metzner (1959), e = ' num2str(case2print.e) ', r = 0 m' ];


% CFD, eccentric case
% Select case2print, last gets indices of cell array elements w/ lengths > 1
case2print = cases(cases.L==0.04 ...
    & cases.e==-0.95 ...    
    & cases.r==3e-5 ...
    & contains(cases.Description,'rpm_p = 0') ...
    & cases.cells==4000 ...
    & strcmp(cases.Fluid,'PAC1') ...
    & strcmp(cases.CurvCorr,'Yes') ...
    & strcmp(cases.IntTrans,''),:);
linetype = '-';
%plot([case2print.U_f_x{:}],-[case2print.dpdx{:}], 'LineStyle',linetype,'LineWidth',2,'Marker','o');

% ScanIT Fig7 of 1st submission
x =    [0.0287	0.336	0.52	0.689	0.893	1.1	1.29	1.47	1.62	1.84]
y = [257	406	673	983	1.43E+03	1.96E+03	2.51E+03	3.06E+03	3.58E+03	4.33E+03]
plot(x,y, 'LineStyle',linetype,'LineWidth',2,'Marker','o');

legendentries{4,1} = [ 'CFD, e = ' num2str(case2print.e) ', r = ' num2str(case2print.r) ' m' ];


% Friction factors
vel = linspace(0.1,2.5)';
Re = case2print.rho_f.*vel.^(2-case2print.n_PL_f).*case2print.d_h.^case2print.n_PL_f./(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f);


% PL, turbulent regime (Irvine 1988)
linetype = '--'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
f=zeros(length(Re),1);
for ii = 1:length(Re)
   if Re(ii)<Re_cr
       f(ii)=24/Re(ii);
   else
       f(ii)=((((2.^(case2print.n_PL_f+4))./(7.^(7.*case2print.n_PL_f))).*(4.*case2print.n_PL_f./(3.*case2print.n_PL_f+1)).^(3.*case2print.n_PL_f.^2))./Re(ii)).^(1./(3.*case2print.n_PL_f+1));
   end
end
f_cr = 24/Re_cr;
f_cr = ((((2.^(case2print.n_PL_f+4))./(7.^(7.*case2print.n_PL_f))).*(4.*case2print.n_PL_f./(3.*case2print.n_PL_f+1)).^(3.*case2print.n_PL_f.^2))./Re_cr).^(1./(3.*case2print.n_PL_f+1));
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
plot(vel,dpdx,'LineStyle',linetype);
legendentries{5,1} = ['Irvine (1988), e = ' num2str(case2print.e) ', r = 0 m' ];


% PL, turbulent regime (Dodge & Metzner 1959), implicit in f hence specify
% existing f >> Re >> vel
linetype = ':'; set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
for ii = 1:length(Re)
   if f(ii)>f_cr
       Re(ii)=24/f(ii);
   else
       Re(ii) = (10.^((1./sqrt(f(ii))+0.4./(case2print.n_PL_f.^1.2)).*(case2print.n_PL_f.^0.75)./4))./(f(ii).^((2-case2print.n_PL_f)./2));
   end
end
vel = (Re.*(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f)./(case2print.rho_f.*case2print.d_h.^case2print.n_PL_f)).^(1./(2-case2print.n_PL_f));
dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
dpdx = EccentricPressureCorrection( dpdx, Re, Re_cr, case2print.d_i, case2print.d_o, -case2print.e, case2print.n_PL_f  );
plot(vel,dpdx,'LineStyle',linetype);
legendentries{6,1} = ['Dodge & Metzner (1959), e = ' num2str(case2print.e) ', r = 0 m' ];



% Experimental data of UiS
set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
expdata = xlsread(filepath_metadata,'Khatibi2018_ExpData','D9:E30');
scatter(expdata(:,1),expdata(:,2));
% scatter(cases2print.U_f_x_exp{aa,:},-cases2print.dpdx_exp{aa,:});
legendentries{7,1} = [ 'Khatibi et al. (2018b), e = ' num2str(cases2print.e(aa)) ', r = ' num2str(case2print.r) ' m' ];



% Dosunmu & Shah (2015)
% It is unclear which Re definition has to be used for the friction factor
% correlation: The MR or the one described in the paper of dosunmu & Shah
% which is the one of Pilevarie % Serth (2009). The to definitions give
% different Re for the same velocity.

% vel = linspace(0.1,2.5)';
% Re = case2print.rho_f.*vel.^(2-case2print.n_PL_f).*case2print.d_h.^case2print.n_PL_f./(12.^(case2print.n_PL_f-1).*case2print.K_PL_f.*((2.*case2print.n_PL_f+1)./(3.*case2print.n_PL_f)).^case2print.n_PL_f);
% 
% kappa = case2print.d_i/case2print.d_o;
% beta = (1 + kappa).*(1-kappa).^(2+1./case2print.n_PL_f).*((1-1./93.*case2print.n_PL_f).^(5./9).*(1./kappa-1).^(9/10)).^(-1);
% d_eff = beta.*case2print.d_o.^((3.*case2print.n_PL_f+1)./case2print.n_PL_f)./((case2print.d_i+case2print.d_o).*(case2print.d_o-case2print.d_i).^((case2print.n_PL_f+1)./case2print.n_PL_f));
% Re = case2print.rho_f.*vel.^(2-case2print.n_PL_f).*d_eff.^case2print.n_PL_f./(12.^(case2print.n_PL_f-1).*case2print.K_PL_f);
% f = 0.002527.*kappa+16.966./Re.^1.1106+0.04102.*case2print.r./case2print.d_h-0.000521;
% dpdx = f.*case2print.rho_f.*vel.^2./(case2print.d_h./2);
% plot(vel,dpdx,'LineStyle',linetype,'Marker','o');


















% Adjust axis
set(gca,...
        'FontSize',12,...
    'xlim', [0 2.5],...
    'ylim', [0 4500]);

% Legend
TightFigure(1);
legendentries = legendentries(~cellfun('isempty',legendentries));
legend(legendentries,'Location','northwest');

% Print to files
fig_name = 'SP_dpdx-U_PAC1';
set(fig_handle(8),'PaperPositionMode','auto','PaperType','A5','PaperOrientation','Landscape') %set paper pos for printing
print(fig_handle(8),[fig_path fig_name],'-dpdf','-bestfit');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(8),[fig_path fig_name],'-depsc');  % '-dpng' '-depsc' 'jpeg' '-depsc2' '-depsc'
print(fig_handle(8),[fig_path fig_name],'-dmeta');  % '-dpng' '-depsc' 'jpeg'


