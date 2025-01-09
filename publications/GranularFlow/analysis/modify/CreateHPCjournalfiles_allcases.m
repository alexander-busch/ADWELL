clc;
close all;
clear all;

% Path to generic journal file
path2journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\hpc';

currentpath = pwd;
cd(path2journal);

% Aspect ratio
AR = [1 2 3]; 

% Type of simulation
% sim = 'run2init';
sim = 'run';

% Fluid
fluid = {'air' 'h2o' 'pac2' 'pac4'};
crosscoefficients  = [0.0109 0.414 0.0721 0.00102;...
    0.0261 0.392 0.21 0.00102];    
    % Cross coefficients
    % Columns K n mu_0 mu_inf
    % Rows pac2 pac4PAC4

% Spatial scale "Domain size"
x_min1 = [-0.05 -0.5 -5];
x_min2 = [-0.04 -0.4 -4];
x_line = [-0.045 -0.45 -4.5];
y_max1 = [0.04 0.4 4];
% y_max_2 is dpendent on AR and is therefore set in the AR loop
alpha_ini = 0.55; %  [0.55 0.575 0.60]

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
    '(rpsetvar ''y_max2 0.01)',...
    '(rpsetvar ''cross/k 1)',...
    '(rpsetvar ''cross/n 1)',...
    '(rpsetvar ''cross/mu_0 0)',...
    '(rpsetvar ''cross/mu_inf 0.001002)',...    
    '(rpsetvar ''sand/alpha_ini 0.55)',...
    '(rpsetvar ''sand/d 0.0001)',...
    'patch_1.ip',...
    'fluid y air',...
    '/solve/set/time-step 1e-4',...
    '/solve/dual-time-iterate 40000 40',...
    ';/define/user-defined/',...
    '; ---- NEW CASE STARTS HERE',...
    'solutions/a=1/alpha_s=0.55/air/1',...
    '"rheology.c"'};

% Number of cases in job array
n = 9;

% Purge
if exist('journals_hpc','dir')
    % rmdir([pwd '\journals_hpc'],'s');
end
mkdir(pwd,'journals_hpc');

% read generic journalfile
fid = fopen([sim '.jou'],'rt') ;
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;

% AR = AR(1); 
% fluid = fluid{2};

for aa = 1:length(AR)
    
    if strcmp(sim,'run2init')
        y_max2 = 1.1111*AR(aa).*(x_min2-x_min1);
    else 
        y_max2 = AR(aa).*(x_min2-x_min1);
    end
    
    
    
    for bb = 1:length(alpha_ini)
        for cc = 1:length(fluid)
            % cc = 2;
            
            % Loop all case combinations and generate one journalfile per case
            index = 1; %Case number
            for ii = 1:length(x_min1)
                for jj = 1:length(d_s)

                    for hh = 1:length(stringlist)
                        % hh=20
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
                                if strcmp(fluid{cc},'pac2')
                                    S2 = ['(rpsetvar ''cross/k ' num2str(crosscoefficients(1,1)) ')'];
                                elseif strcmp(fluid{cc},'pac4')
                                    S2 = ['(rpsetvar ''cross/k ' num2str(crosscoefficients(2,1)) ')'];
                                else
                                    S2 = ';No cross model coefficients required.';
                                end
                            case 10
                                if strcmp(fluid{cc},'pac2')
                                    S2 = ['(rpsetvar ''cross/n ' num2str(crosscoefficients(1,2)) ')'];
                                elseif strcmp(fluid{cc},'pac4')
                                    S2 = ['(rpsetvar ''cross/n ' num2str(crosscoefficients(2,2)) ')'];
                                else
                                    S2 = ';No cross model coefficients required.';
                                end
                            case 11
                                if strcmp(fluid{cc},'pac2')
                                    S2 = ['(rpsetvar ''cross/mu_0 ' num2str(crosscoefficients(1,3)) ')'];
                                elseif strcmp(fluid{cc},'pac4')
                                    S2 = ['(rpsetvar ''cross/mu_0 ' num2str(crosscoefficients(2,3)) ')'];
                                else
                                    S2 = ';No cross model coefficients required.';
                                end
                            case 12
                                if strcmp(fluid{cc},'pac2')
                                    S2 = ['(rpsetvar ''cross/mu_inf ' num2str(crosscoefficients(1,4)) ')'];
                                elseif strcmp(fluid{cc},'pac4')
                                    S2 = ['(rpsetvar ''cross/mu_inf ' num2str(crosscoefficients(2,4)) ')'];
                                else
                                    S2 = ';No cross model coefficients required.';
                                end
                            case 13
                                S2 = ['(rpsetvar ''sand/alpha_ini ' num2str(alpha_ini(bb)) ')'];
                            case 14
                                S2 = ['(rpsetvar ''sand/d ' num2str(d_s(jj)) ')'];
                            case 15
                                S2 = ['patch_' num2str(index) '.ip'];
                            case 16
                                if (strcmp(fluid{cc},'pac2') || strcmp(fluid{cc},'pac4'))
                                    S2 = 'fluid y cross';
                                else
                                    S2 = ['fluid y ' fluid{cc}];                 
                                end
                            case 17
                                if strcmp(fluid{cc},'air')
                                    S2 = '/solve/set/time-step 1e-4';
                                else
                                    S2 = '/solve/set/time-step 1e-3';
                                end
                            case 18
                                if strcmp(fluid{cc},'air')
                                    S2 = '/solve/dual-time-iterate 20000 40';
                                elseif strcmp(fluid{cc},'h2o')
                                    if (alpha_ini==0.55)
                                        S2 = '/solve/dual-time-iterate 20000 40';
                                    else
                                        S2 = '/solve/dual-time-iterate 60000 40';
                                    end
                                elseif strcmp(fluid{cc},'pac2')
                                    if (alpha_ini==0.55)
                                        S2 = '/solve/dual-time-iterate 40000 40';
                                    else
                                        S2 = '/solve/dual-time-iterate 120000 40';
                                    end
                                else
                                    if (alpha_ini==0.55)
                                        S2 = '/solve/dual-time-iterate 60000 40';
                                    else
                                        S2 = '/solve/dual-time-iterate 130000 40';
                                    end
                                end
                            case 19
                                if (strcmp(fluid{cc},'pac2') || strcmp(fluid{cc},'pac4'))
                                    S2 = '/define/user-defined/';
                                else
                                    S2 = ';No UDFs required ...';
                                end
                            case 20
                                S2 = ['; ---- NEW CASE STARTS HERE: a = ' num2str(AR(aa)) ', alpha_s = ' num2str(alpha_ini(bb)) ', ' fluid{cc} ', case' num2str(index)];
                            case 21
                                S2 = ['solutions/a=' num2str(AR(aa)) '/alpha_s=' num2str(alpha_ini(bb)) '/' fluid{cc} '/' num2str(index)];
                            case 22
                                if strcmp(fluid{cc},'pac2')
                                    S2 = '"rheology_pac2.c"';
                                elseif strcmp(fluid{cc},'pac4')
                                    S2 = '"rheology_pac4.c"';
                                else
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

                    end % of replace strings loop
                    

                    % Assemble journal file content
                    if index==1
                        Z = Y;
                    else
                        Z = [Z Y];
                    end

                    index=index+1;
                end % of d_s loop
                
            end % of x_min1 loop
            
        
            % Assemble final journal file content
            if ( cc==1 )
                Z_cc = Z;
            else
                Z_cc = [Z_cc Z];
            end   
        
        end % of fluid loop
        
        
        % Assemble final journal file content
        if ( (bb==1) )
            Z_bb = Z_cc;
        else
            Z_bb = [Z_bb Z_cc];
        end
        
    end % of alpha_ini loop
    
       
    % Assemble final journal file content
    if aa==1
        Z_aa = Z_bb;
    else
        Z_aa = [Z_aa Z_bb];
    end
        
end % of AR loop


% Remove recorded journal file
fid_remove = fopen('journals_archive/yaxislog.jou','rt');
X_remove = fread(fid_remove);
fclose(fid_remove);
X_remove = char(X_remove.');
Z_aa = strrep(Z_aa, X_remove, '; Removed because incompatible with hpc');

% 
Z_aa = strrep(Z_aa, 'Paste "/file/read-journal run.jou"', 'Paste "/file/read-journal run_all.jou"');


% For testing of journal file
% Z_aa = strrep(Z_aa, '/solve/dual-time-iterate 40000 40', '/solve/dual-time-iterate 40 40');
% Z_aa = strrep(Z_aa, '/solve/dual-time-iterate 160000 40', '/solve/dual-time-iterate 100 40');


% Create new journalfile
% fid2 = fopen('journals_hpc/run_all.jou','wt') ;
fid2 = fopen('run_all.jou','wt') ;
fwrite(fid2,Z_aa) ;
fclose (fid2) ;


cd(currentpath);
