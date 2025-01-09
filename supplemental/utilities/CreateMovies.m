% Read png images from path and creates mp4 movies


%% Cleanup
close all
clear all
clc

%% Specification

% CFD cases
casenamelist = {'161124_vp_pac1037_s100',...
    '161124_vp_pac1037_s200',...
    '161124_vp_pac2037_s100',...
    '161124_vp_pac2037_s200',...
    '161124_vp_pac3037_s200'};

% List of different image types
imagenamelist = {'xy_fluid_vel',...
    'xy_solid_vel',...
    'xy_solid_vof',...    
    'yz_fluid_vel',...
    'yz_solid_vel',...
    'yz_solid_vof',...
    'xyz_fluid_vel',...
    'xyz_solid_vel',...
    'xyz_solid_vof'};

%     'front_fluid_vel',...
%     'front_solid_vel',...
%     'front_solid_vof'};

%% Create videos
for i= 1:length(casenamelist) % Loop through all CFD cases

% Assemble path to images
pathname = ['C:\Users\alexabus\IEPT1325_Local\Data\SW_working_directories\FLUENT\HPC\Simulations\' casenamelist{i} '\solution'];

% Get amount of images per type
% images = length(dir([pathname '\images\' imagenamelist{1} '*.png']));   

    for j = 1:length(imagenamelist) % Loop through various image types
        
        % Delete movie if exists
        moviepath = [pathname '\movies\' imagenamelist{j} '.mp4']; % Create path to movie
        if exist(moviepath, 'file')==2
            delete(moviepath);
        end
        
        % Create movie
        list = dir([pathname '\images\' imagenamelist{j} '*.png']); % Create structure off all files where the first column is the filename
        writerObj = VideoWriter([pathname '\movies\' imagenamelist{j} '.mp4'],'MPEG-4'); % Create video object
        writerObj.FrameRate = 10; % Specify frame rate
        open(writerObj); % Open video object

        for k = 1:numel(list); % Loop through all images of type imagenamelist{i};
            filename = list(k).name; % Create filename of current image
            thisimage = imread([pathname '\images\' filename],'png'); % Read current image 
            writeVideo(writerObj, thisimage); % Write current image to video object
        end

        close(writerObj); % Close video object
        
        disp ([imagenamelist{j} '.mp4']);
    end

end
