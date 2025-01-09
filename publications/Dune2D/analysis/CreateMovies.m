% List of different image types
imagenamelist = {'xy_fluid_vel',...
    'xy_solid_vel',...
    'xy_solid_vof'};

% Plot name list
plotnamelist = {'Fluid velocity',...
    'Solid velocity',...
    'Solid volume fraction'};

% Pictures per RT second
framerate = 20; 

text_str = cell(5,1);
text_str{1} = 'Water & sand 2D channel flow, ';
text_str{2} = 'Solid density 2560 kg/m³, Particle diameter 0.0003 m, Fluid superficial velocity 0.43 m/s, Solid superficial velocity 0.43 m/s, Solid volume fraction at inlet 0.015';
text_str{3} = 'Initial bed height 0.005 m, Channel height 0.04 m, Channel length 1 m';
text_str{4} = 'Fluent R17.2, Eulerian-Eulerian 2 fluid, Newtonian, k-omega SST, Segregated, pbns, dx = 0.005 m, dy = 0.002 m, dt = 0.01 s';
text_str{5} = 'Dune2D_170130_1 | Alexander Busch, NTNU | 10.02.2017';

line_pos = [41 74 1560 74; 41 135 1560 135];
text_pos = [36 20; 36 40; 36 140; 36 160; 1300 170];
text_size = [15; 13; 13; 13; 10];

%     thisimage = imread('test.png'); % Read current image  
%     thisimage = insertShape(thisimage,'line',line_pos,'LineWidth',1,'Color','black');
%     for l = 1:length(text_str)
%         thisimage = insertText(thisimage, text_pos(l,:), text_str{l},'FontSize',text_size(l),'BoxColor','white','BoxOpacity',0,'TextColor','black');
%     end
%     imshow(thisimage);

for j = 1:length(imagenamelist) % Loop through various image types

pathname = ['C:\Users\alexabus\IEPT1325_Local\Data\SWwd\FLUENT\AdWell\Dune2D\Model\170130_1 - MP + P + Corrected Source Terms, reworked code\images\'];

list = dir([pathname imagenamelist{j} '*.png']); % Create structure off all files where the first column is the filename
writerObj = VideoWriter([pathname imagenamelist{j} '.mp4'],'MPEG-4'); % Create video object
writerObj.FrameRate = framerate; % Specify frame rate
open(writerObj); % Open video object

for k = 1:numel(list); % Loop through all images of type imagenamelist{i};
	filename = list(k).name; % Create filename of current image
    thisimage = imread([pathname filename],'png'); % Read current image
    thisimage = insertShape(thisimage,'line',line_pos,'LineWidth',1,'Color','black');
    for l = 1:length(text_str)
        if l==1
            text = [text_str{l} plotnamelist{j}];
        else
            text = text_str{l} ;
        end 
        thisimage = insertText(thisimage, text_pos(l,:), text,'FontSize',text_size(l),'BoxColor','white','BoxOpacity',0,'TextColor','black');
    end
    
    writeVideo(writerObj, thisimage); % Write current image to video object
end

close(writerObj); % Close video object
end