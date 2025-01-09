close all;
clear all;
clc;

% Path to folder with images
imagepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\hpc\drillpipemotion';

filename = [imagepath '\drillpipemotion.gif'];

% Get structure of all files for current figure in figurelist
files = dir([imagepath '\*.png']);
    
for ii = 1:length(files)

    % Read image        
    pic = imread([imagepath '\' files(ii).name]);
    
    [imind,cm] = rgb2ind(pic,256);
    
    % Write to the GIF File 
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end   
end