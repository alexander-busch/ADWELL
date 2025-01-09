clc;
close all;
clear all;

% Path to generic journal file
path2journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\hpc';

currentpath = pwd;
cd(path2journal);


% Drill pipe rotation
rpm = [0 30 60 100 130];

% Whirl
whirl = 'yes';

% List of strings to modify
stringlist = {'(rpsetvar ''rpm_pipe 30)'};

% Purge
if exist('journals_hpc','dir')
    rmdir([pwd '\journals_hpc'],'s');
end
mkdir(pwd,'journals_hpc');

% read generic journalfile
% if strcmp(whirl, 'no')
%     fid = fopen('run.jou','rt') ;
% else
%     %fid = fopen('run_and_replace_mesh_serial.jou','rt') ;
%     fid = fopen('run_and_replace_mesh_parallel.jou','rt') ;
% end

fid = fopen('run.jou','rt') ;

X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;


for ii = 1:length(rpm)
    
    for hh = 1:length(stringlist)
        % hh=20
        % String to replace
        S1 = stringlist{hh};
        % New string
        switch hh
            case 1
                S2 = ['(rpsetvar ''rpm_pipe ' num2str(rpm(ii)) ')'];
            case 2

            otherwise
                disp('Error');   
        end

    end % of replace strings loop
    
    % Replace strings
    if hh==1
        Y = strrep(X, S1, S2);
    else
        Y = strrep(Y, S1, S2);
    end  
            
    % Create new journalfile
    fid2 = fopen(['journals_hpc/run' num2str(ii) '.jou'],'wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;   
    
end

cd(currentpath);
