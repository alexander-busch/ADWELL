clear all
clc

HPCpath = 'C:\Users\alexabus\IEPT1325_Local\Data\SW_working_directories\FLUENT\HPC\Simulations\';

foldernamelist = {'vp_h2o1037_s050',...
    'vp_h2o1037_s100',...
    'vp_h2o1037_s200',...
    'vp_pac1037_s100',...
    'vp_pac1037_s200',...
    'vp_pac2037_s100',...
    'vp_pac2037_s200',...
    'vp_pac3037_s200'};

currentFolder = pwd;
cd (HPCpath);

%% Create folders
for i= 1:length(foldernamelist) % Loop through all CFD cases

mkdir (foldernamelist{i})

cd (foldernamelist{i})

mkdir pre
mkdir solution
mkdir temp
 
cd pre
mkdir macros
mkdir reports
mkdir UDFs
 
cd ..
cd solution
mkdir images
mkdir monitors
mkdir movies
mkdir plots
mkdir reports
cd ..
cd ..

end


cd (currentFolder);