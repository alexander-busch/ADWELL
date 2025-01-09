clc;
close all;
clear all;

% Paths
currentpath = pwd;
parentpath = cd(cd('..'));
archivepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive';


% Identifier and fluid (See Metadata excel file for overview of cases)
ident = '181213_1';

% fluid = 'Air (1.225 kg/m³)';
% fluid = 'H2O (1000 kg/m³)';
 fluid = 'PAC2 (1000 kg/m³)';
%fluid = 'PAC4 (1000 kg/m³)';

% Initial solid volume fraction
 alpha_s0 = 0.55;
% alpha_s0 = 0.60;

% Aspect ratio
%  AR = 1;
%AR = 2;
 AR = 3;
%
% Spatial scale "Domain size"
x_max = [0.1 1 10];
y_max = [0.04 0.4 4];
dx = [0.002 0.02 0.2];

% Spatial scale "Particle size"
d_s = [0.0001 0.001 0.01];


% List of strings to modify
stringlist = {'180604_1/1/solution/cfd-post/case-',...
    'case1_t',...
    'sandpile_small.cst',...
    'Text String = Case 1 (W = 0.1 m, H = 0.04 m, dx = 0.002 m, d_s = 0.0001 m)'...
    'Text String = Sand (2650 kg/m³) in air (1.225 kg/m³)',...
    '/180604_1/1/solution/180604_1-1',...
    'Aspect ratio a = 1',...
    'case-100.000000.cdat'};

% read generic file
fid = fopen([archivepath '\cliffcollapse_hpc_createVideo_generic.cse'],'rt');
X = fread(fid);
fclose(fid);
X_generic = char(X.');

% file_beginning = extractBefore(X_generic,'>load filename');
% file_end = extractAfter(X_generic,'# starts here ---------------------------------');



index = 1;
for ii = 1:length(x_max)
    for jj=1:length(d_s)
%         if index==1
%             Y_final = X_generic;
%         else
% 
%         end

        for hh = 1:length(stringlist)
            % String to replace
            S1 = stringlist{hh};
            % New string
            switch hh
                % Debugging hh = 5
                case 1
                    % S2 = [ident '/' num2str(index) '/solution/cfd-post/case-'];
                    S2 = [ident '/' num2str(index) '/cfd-post/case-'];
                case 2
                    S2 = ['case' num2str(index) '_t'];
                case 3
                    if ii==1
                        S2 = 'sandpile_small.cst';
                    elseif ii==2
                        S2 = 'sandpile_medium.cst';
                    else
                        S2 = 'sandpile_large.cst';
                    end                        
                case 4
                    S2 = ['Text String = Case ' num2str(index) ' (W = ' num2str(x_max(ii)) ' m, H = ' num2str(y_max(ii)) ' m, dx = ' num2str(dx(ii)) ' m, d_s = ' num2str(d_s(jj)) ' m)'];
                case 5
                    S2 = ['Text String = Sand (2650 kg/m³) in ' fluid];
                case 6
                    S2 = ['/' ident '/' ident '-' num2str(index)];
                case 7
                    S2 = ['Aspect ratio a = ' num2str(AR) ', alpha_s,0 = ' num2str(alpha_s0) ];
                case 8
                    if (alpha_s0 == 0.55)
                        if contains(fluid,'Air')
                            S2 = 'case-2.000000.cdat';
                        elseif contains(fluid,'H2O')
                            S2 = 'case-20.000000.cdat';
                        elseif contains(fluid,'PAC2')
                            S2 = 'case-40.000000.cdat';
                        else
                            S2 = 'case-60.000000.cdat';
                        end
                    else
                        if contains(fluid,'Air')
                            S2 = 'case-2.000000.cdat';
                        elseif contains(fluid,'H2O')
                            S2 = 'case-60.000000.cdat';
                        elseif contains(fluid,'PAC2')
                            S2 = 'case-120.000000.cdat';
                        else
                            S2 = 'case-130.000000.cdat';
                        end
                    end
                otherwise
                    disp('Error');
            end
            % Replace strings
            if hh==1
                Y = strrep(X_generic, S1, S2);
            else
                Y = strrep(Y, S1, S2);
            end
        end
        
        
        
        if index==1
            Y_final = Y;
        else
            Y_final = [Y_final Y];
        end

        index=index+1;
    end
end




%% Remove 

fid_remove = fopen([pwd '\remove_CreateVideos.txt'],'rt');
X = fread(fid_remove);
fclose(fid_remove);
X_remove = char(X.');
Y_final = strrep(Y_final, X_remove, '');

% Create new sessionfile
fid2 = fopen([archivepath '\cliffcollapse_hpc_createVideo.cse'],'wt');
fwrite(fid2,Y_final);
fclose (fid2);

% cd(path2journal);
% cd(currentpath);

