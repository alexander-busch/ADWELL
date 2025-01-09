% Change folder names


dirName = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive\180425_1';

dirResult = dir(dirName);
allDirs = dirResult([dirResult.isdir]);
allSubDirs = allDirs(3:end);

for ii = length(allSubDirs):-1:2
    thisDir = allSubDirs(ii);
    oldname = fullfile(dirName,thisDir.name);
    newname = [dirName '\' num2str(ii)];
        movefile(oldname,newname);
end

% for i = 1:length(allSubDirs)
%     thisDir = allSubDirs(i);
%     thisDirName = thisDir.name;
%     if ~strcmp(thisDirName(end-2:end),'1-4')
%         oldname = fullfile(dirName,thisDir.name);
%         newname = [fullfile(dirName,thisDir.name) '1-4'];
%         movefile(oldname,newname);
%     end
%     changeDirNames(newname);
% end
% 



