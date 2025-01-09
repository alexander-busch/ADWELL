clc;
close all;
clear all;

% Identifier and fluid (See Metadata excel file for overview of cases)
ident = '190403_1';
fluid = 'pac4'; % air h2o pac2 pac4 
a = 3; % 1 2 3
alpha_s0 = 0.60; % 0.55 0.60

% Paths
currentpath = pwd;
parentpath = cd(cd('..'));
archivepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive';

% Purge and create folder for CFD post exports
folder = [parentpath '\cfdpostexports\alpha_s=' num2str(alpha_s0) '\a=' num2str(a) '\' fluid];
if exist(folder,'dir')
    rmdir(folder,'s');
end
mkdir(folder);

% Spatial scale "Domain size"
x_max = [0.1 1 10];
y_max = [0.04 0.4 4];
dx = [0.002 0.02 0.2];

% Spatial scale "Particle size"
d_s = [0.0001 0.001 0.01];

% List of strings to modify
stringlist = {'180604_1/1/solution/',...
    'cfdpostexports/case1_t',...
    'sandpile_small.cst'};


% read generic file
if (alpha_s0 == 0.55)
    if(strcmp(fluid,'air') || strcmp(fluid,'h2o'))
        % fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_air.cse'],'rt');
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_air2.cse'],'rt');
        % Air has dt = 1e-4 s and t = 0...2s and h2o has dt = 1e-3 and t =
        % 0...20s, therefore both use the ...air2.cse file where the 10000th and
        % 20000th are read into CFD-post.
    elseif strcmp(fluid,'pac2')
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_0.55_pac2.cse'],'rt');
    else
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_0.55_pac4.cse'],'rt');
    end
    
else
    if strcmp(fluid,'air')
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_air2.cse'],'rt');
    elseif strcmp(fluid,'h2o')
        % Same flowtime as PAC 4 0.55 
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_0.55_pac4.cse'],'rt');
    elseif strcmp(fluid,'pac2')
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_0.60_pac2.cse'],'rt');
    else
        fid = fopen([archivepath '\cliffcollapse_hpc_exportVOF_generic_0.60_pac4.cse'],'rt');
    end
    
end

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
                case 1
                    % S2 = [ident '/' num2str(index) '/solution/'];
                    S2 = [ident '/' num2str(index) '/'];
                case 2
                    % S2 = ['cfdpostexports/a=' num2str(a) '/' fluid '/case' num2str(index) '_t'] % ,'%02d' provides leading zeros
                    S2 = ['cfdpostexports/alpha_s=' num2str(alpha_s0) '/a=' num2str(a) '/' fluid '/case' num2str(index) '_t'] % ,'%02d' provides leading zeros
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

fid_remove = fopen([pwd '\remove_ExportVOF.txt'],'rt');
X = fread(fid_remove);
fclose(fid_remove);
X_remove = char(X.');
Y_final = strrep(Y_final, X_remove, '');

% Create new sessionfile
fid2 = fopen([archivepath '\cliffcollapse_hpc_exportVOF.cse'],'wt');
fwrite(fid2,Y_final);
fclose (fid2);

% cd(path2journal);
% cd(currentpath);

