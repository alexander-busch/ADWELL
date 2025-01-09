close all;
clear all;
clc;


%addpath(genpath('C:\Users\U376949\OneDrive - Danfoss\Documents\MATLAB'));

% Path to folder with images
imagepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive\181027_2\5\solution\images';


cases = {'181026_1',...
    '181026_2',...
    '181026_3',...
    '181027_1',...
    '181027_2',...
    '181121_1',...
    '181121_2',...
    '181121_3',...
    '181028_1',...
    '181028_2'};

mkdir([imagepath(1:end-6) 'movies']);

% List of various figures 
figurelist = {'fluid-velocities_solid-volume-fraction_xy';...
    'fluid-velocities_solid-volume-fraction_yz';...
    'fluid-velocities_solid-volume-fraction_iso';...
    'solid-velocities_solid-volume-fraction_xy';...
    'solid-velocities_solid-volume-fraction_yz';...
    'solid-velocities_solid-volume-fraction_iso'};

% Pressure gradient profile
dpdx = (50:100:550)';
dpdx = (0:100:500)';
dt = 1e-3;
tmax = GetRealFlowTime( 30 );
tmax = 12*2.374; % Fluent
idx_up = tmax*(1/dt);
time = dt*(0:1:length(dpdx)-1)'*idx_up;

% Loop cases
for aa=1:length(cases)
    imagepath = ['C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\archive\' cases{aa} '\5\solution\images'];
    moviepath = [imagepath(1:end-6) 'movies'];
    if exist(moviepath,'dir')
        rmdir(moviepath)
    end
    
    mkdir(moviepath);


% Loop filelist
    for ii=1:length(figurelist)
        % ii=3
        % Create video object
        video = VideoWriter([moviepath '\' figurelist{ii} ],'MPEG-4');

        % Frame rate
        video.FrameRate = 33;

        % Open video object
        open(video);

        % Get structure of all files for current figure in figurelist
        files = dir([imagepath '\' figurelist{ii} '_' '*.png']);

        % Loop structure
        for jj=1:length(files)

            % Read image        
            pic = imread([imagepath '\' files(jj).name]);

            % Convert to double
            pic = im2double(pic);

            % Get time
            currenttime = round(str2double(files(jj).name(end-13:end-4)),3);

            % Get dpdx as function of time

            currentdpdx = interp1(time,dpdx,currenttime,'next');

            % Get image size, second value is x and first value is y,
            % coordinate system starts at top left corner
            pic_size = size(pic);  

            % Add dpdx and time
            pic = insertText(pic,[pic_size(2) 0],['t = ' files(jj).name(end-13:end-4) ' s'],...
                'FontSize',18,...
                'AnchorPoint','RightTop',...
                'BoxColor','w');


            pic = insertText(pic,[pic_size(2) 50],['dp/dx = ' num2str(currentdpdx) ' Pa/m'],...
                'FontSize',18,...
                'AnchorPoint','RightTop',...
                'BoxColor','w');

            % imshow(pic)

            % Write to video object
            writeVideo(video,pic);
        end

        % Close video object
        close(video);

    end
end