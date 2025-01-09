clear all
clc


casenamelist = {'161124_vp_pac1037_s100',...
    '161124_vp_pac1037_s200',...
    '161124_vp_pac2037_s100',...
    '161124_vp_pac2037_s200',...
    '161124_vp_pac3037_s200'};


for i = 1:length(casenamelist)

    % Assemble path to images
    d = ['C:\Users\alexabus\IEPT1325_Local\Data\SW_working_directories\FLUENT\HPC\Simulations\' casenamelist{i} '\solution\images'];

    names = dir(d);
    names = {names(~[names.isdir]).name};

    for n = 1:numel(names)
        %if length(names{n})<24 %including extension
        oldname = [d '\' names{n}]; % Old filename incl path

        a=strfind(oldname,'-'); % Get position of - in string
        b=strfind(oldname,'.'); % Get position of . in string

        if b(1)-a == 2
            if strfind(oldname,'_vof-')>0
                C = strsplit(oldname,'_vof-');    
                newname = strcat(C(1),'_vof-0',C(2));
            else
                C = strsplit(oldname,'_vel-');    
                newname = strcat(C(1),'_vel-0',C(2));
            end
            movefile(oldname,newname{1});
        end
    end
end

% Alternative to movefile

 
%doscomand = strcat('rename', {' '}, oldname, {' '}, newname);
%dos(doscomand{1}); % (1)

% newname must not include the full path, just the newname.
