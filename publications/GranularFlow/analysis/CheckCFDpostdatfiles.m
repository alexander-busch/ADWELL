clc;
close all;
clear all;

addpath 'C:\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange';

fluid = 'water'; % water pac2 pac4

datapath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive\180604_1';

for aa = 1:9
    for bb = 1:10
        bbb=10*bb;
        if exist([ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb) '.000000.cdat'],'file')
            disp(['Case ' num2str(aa) ', file ' num2str(bbb) '.000000.cdat ok']);
        else
            if exist([ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb) '.001000.cdat'],'file')
                
                copyfile([ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb) '.001000.cdat'],...
                    [ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb) '.000000.cdat']);
                disp(['Case ' num2str(aa) ', file ' num2str(bbb) '.000000.cdat created from file ' num2str(bbb) '.001000.cdat']);
                
            elseif exist([ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb-1) '.999000.cdat'],'file')
                
                copyfile([ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb-1) '.999000.cdat'],...
                    [ datapath '\' num2str(aa) '\solution\cfd-post\case-' num2str(bbb) '.000000.cdat']);
                disp(['Case ' num2str(aa) ', file ' num2str(bbb) '.000000.cdat created from file ' num2str(bbb-1) '.999000.cdat']);
                
            else
                
                disp(['Case ' num2str(aa) ', file ' num2str(bbb) '.000000.cdat ERROR']);
                
            end
                
        end
    end
end