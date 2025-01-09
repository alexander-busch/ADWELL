clc;
close all;
clear all;

% Path to generic journal file
path2journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\hpc';

currentpath = pwd;
cd(path2journal);

% Aspect ratio
a = 2;
% a = 3;

% Type of simulation
% sim = 'run2init';
sim = 'run';

% Fluid
fluid = 'air'; % air h2o cross

% Spatial scale "Domain size"
x_min1 = [-0.05 -0.5 -5];
x_min2 = [-0.04 -0.4 -4];
x_line = [-0.045 -0.45 -4.5];
y_max1 = [0.04 0.4 4];
if strcmp(sim,'run2init')
    y_max2 = 1.1111*a.*(x_min2-x_min1);
else 
    y_max2 = a.*(x_min2-x_min1);
end
alpha_ini = [0.55 0.55 0.55];

% Spatial scale "Particle size"
d_s = [0.0001 0.001 0.01];

% List of strings to modify
stringlist = {'case_1_2patch.cas',...
    'case_1.cas',...
    '(rpsetvar ''x_min1 -0.05)',...
    '(rpsetvar ''x_min2 -0.04)',...
    '(rpsetvar ''x_line -0.045)',...
    '(rpsetvar ''x_max 0.05)',...
    '(rpsetvar ''y_max1 0.04)',...
    ['(rpsetvar ''y_max2 ' num2str(y_max2(1)) ')'],...
    '(rpsetvar ''sand/alpha_ini 0.55)',...
    '(rpsetvar ''sand/d 0.0001)',...
    'patch_1.ip',...
    'fluid y air',...
    '/solve/set/time-step 1e-4',...
    '/solve/dual-time-iterate 40000 40'};

% Number of cases in job array
n = 9;

% Purge
if exist('journals_hpc','dir')
    rmdir([pwd '\journals_hpc'],'s');
end
mkdir(pwd,'journals_hpc');

% read generic journalfile
fid = fopen([sim '.jou'],'rt') ;
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;


% Loop all case combinations and generate one journalfile per case
index = 1; %Case number
for ii = 1:length(x_min1)
    for jj = 1:length(d_s)

        for hh = 1:length(stringlist)
            % String to replace
            S1 = stringlist{hh};
            % New string
            switch hh
                case 1
                    S2 = ['case_' num2str(ii) '_2patch.cas'];
                case 2
                    S2 = ['case_' num2str(ii) '.cas'];
                case 3
                    S2 = ['(rpsetvar ''x_min1 ' num2str(x_min1(ii)) ')'];
                case 4
                    S2 = ['(rpsetvar ''x_min2 ' num2str(x_min2(ii)) ')'];
                case 5
                    S2 = ['(rpsetvar ''x_line ' num2str(x_line(ii)) ')'];
                case 6
                    S2 = ['(rpsetvar ''x_max ' num2str(-x_min1(ii)) ')'];
                case 7
                    S2 = ['(rpsetvar ''y_max1 ' num2str(y_max1(ii)) ')'];
                case 8
                    S2 = ['(rpsetvar ''y_max2 ' num2str(y_max2(ii)) ')'];
                case 9
                    S2 = ['(rpsetvar ''sand/alpha_ini ' num2str(alpha_ini(ii)) ')'];
                case 10
                    S2 = ['(rpsetvar ''sand/d ' num2str(d_s(jj)) ')'];
                case 11
                    S2 = ['patch_' num2str(index) '.ip'];
                case 12
                    S2 = ['fluid y ' fluid];
                case 13
                    if strcmp(fluid,'air')
                        S2 = '/solve/set/time-step 1e-4';
                    else
                        S2 = '/solve/set/time-step 1e-3';
                    end
                case 14
                    if strcmp(fluid,'air')
                        S2 = '/solve/dual-time-iterate 40000 40';
                    else
                        S2 = '/solve/dual-time-iterate 160000 40';
                    end
                otherwise
                    disp('Error');
                    
                    
                    
                    
            end
            % Replace strings
            if hh==1
                Y = strrep(X, S1, S2);
            else
                Y = strrep(Y, S1, S2);
            end          
        end
        
        % Remove recorded journal file
        fid_remove = fopen(['journals_archive/yaxislog.jou'],'rt');
        X_remove = fread(fid_remove);
        fclose(fid_remove);
        X_remove = char(X_remove.');
        Y = strrep(Y, X_remove, '; Removed because incompatible with hpc');

        % Create new journalfile
        fid2 = fopen(['journals_hpc/run' num2str(index) '.jou'],'wt') ;
        fwrite(fid2,Y) ;
        fclose (fid2) ;

        index=index+1;
    end
   
    % movefile(['journalfiles/run' num2str(i) '.jou'],[pwd '/run' num2str(i) '.jou']);
end

cd(currentpath);
