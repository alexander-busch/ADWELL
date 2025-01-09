% Re analysis

clear all;
close all;
clc;


Import;
Clean;

% addpath('M:\Documents\SW_working_directories\MATLAB\ParaSpace');
% load ParaSpace.mat;


%% Plot

figure;
hold on;

% Format
xlabel('Hole [in]');
ylabel('Re [-]');
grid('on');
set(gca,...
    'YScale','log',...
    'box','on',...
    'FontSize',24);
set(gcf,...
    'color','w');



Re_Range = linspace(0,200000,2000); % Definition of Re range



for i=1:length(D_H)
Re = kron(Re_pre,D_H(i));
% plot(D_H(i),min(min(Re)),'ko');
% h = histogram(Re,RE_Range,'Visible','Off');
% bar3(h.Values,D_H(i));
    for j=1:size(Re,1)
        for k=1:size(Re,2)
            scatter(Hole(i),Re(j,k),'ko');
        end
    end
end



%
