% Parameters



%% Geometry [in]
% Typical combinations of diameter values for annulus (d_o and d_i) and
% drill pipe only (d_i = 0)[in]
% See https://www.woolleytool.com/files/recommended_drill_pipe.pdf
% First row is the inner drill pipe diameter,
% second row is the bit/hole diameter and the outer drill pipe diameter
d_o = in2m*[5.965 5.965 5.965 5.965 4.670 4.21401 4.214 3.826;...
    42.000 36.000 26.000 17.500 12.250 09.875 08.500 06.625]';
d_i = in2m*[00.000 00.000 00.000 00.000 00.000 00.000 00.000 00.000;...
    06.625 06.625 06.625 06.625 05.500 05.000 05.000 04.500]';
% Redefinition based on recommendations from TAG
d_o = in2m*[5.965 5.965 5.965 5.965 4.670 4.214 3.826;...
    42.000 36.000 26.000 17.500 12.250 08.500 06.5]';
d_i = in2m*[00.000 00.000 00.000 00.000 00.000 00.000 00.000;...
    06.625 06.625 06.625 06.625 05.500 05.000 04.500]';

% Only inclined sections
% d_o = in2m*[5.965 4.670 4.214 3.826;...
%     17.500 12.250 08.500 06.5]';
% d_i = in2m*[00.000 00.000 00.000 00.000;...
%     06.625 05.500 05.000 04.500]';






sections = categorical(d_o(:,2)/in2m); % Categorical array of Wellbore section for plotting bar charts
% Hydraulic diameter
d_h = d_o-d_i;
% Cross-sectional area
A = pi/4*(d_o.^2-d_i.^2);
% Create cell array containing d_h, d_i, d_o for x-axis
% Matrix with one string per section
matrix=cell(length(d_o)+2,1);
for i=1:length(matrix)
    if (i==1 || i==length(matrix))        
        matrix{i} = '';
    else
        matrix{i} = ['\color{black}',num2str(d_h(i-1,2),2),'\color{black}, '...
            '\color{black}',num2str(d_h(i-1,1),2),...
            '\newline       ','\color{black}',num2str(d_i(i-1,2)/in2m),...
            '\newline       ',num2str(d_o(i-1,2)/in2m)];
    end
end


% Matrix with individual strings per section, issue: onlz last is displayed
% on axis
matrix2=cell(length(d_o)+2,3);
for i=1:length(matrix2)
    if (i==1 || i==length(matrix2))        

    else
        matrix2{i,1} = ['\color{black}',num2str(d_h(i-1,2),2),'\color{black}, '...
            '\color{black}',num2str(d_h(i-1,1),2)];
        matrix2{i,2} = ['\color{black}',num2str(d_i(i-1,2)/in2m)];
        matrix2{i,3} = [num2str(d_o(i-1,2)/in2m)];
    end
end


%% Flow rates [gpm]
% !!! Check for large hole sections and update 2000 and 2500 !!!
% first column is the min. flow rate, second row is the max. flow rate as
% recommended by K&M (Mims, M., Krepp, T., 2007. Drilling  Design  and
% Implementation  For  Extended  Reach  and  Complex Wells, 3rd ed. K&M Technologies Group, LLC.)
Q = gpm2m3ps*[2000 2000 2000 900 750 350 150;...
    2500 2500 2500 1500 1100 600 250]';
%Q = [Q Q]'; % Q(:,4:length(Q))];
% 
% Q = gpm2m3ps*[900 750 350 150;...
%     1500 1100 600 250]';


%% Drill string rotation rates [rpm]
% !!! Check for large and small hole sections and update 300 and 70 !!!
rpm = [240 220 200 180 120 70 70]';
rpm = [rpm 2*rpm];



%% Material properties

% Particle diameter
d_p = [0.1 1 10]/1000;

% Densities
rho_f = 1000;
rho_s = 2650;

% Drilling fluid rheological properties, representative oilfield PV & YP acc. to Busch et al. (2016)
YP = [0; 10; 20]; % [Pa]
PV = [15; 30; 50]/1000; % [Pa.s]
% Assemble vectors for YP and PV
YP = [ ones(length(YP),1)*YP(1); ones(length(YP),1)*YP(2); ones(length(YP),1)*YP(3)];
PV = [ PV; PV; PV]; % [mPa.s]
% Add water
YP = [ 0 YP' ]';
PV = [ 0.00102 PV']';

% Shear rate range
SR_max = 10000;
SR_range = logspace(-2,log10(SR_max));
SR_Fann = [600 300 200 100 6 3]*FannRPM2SR;