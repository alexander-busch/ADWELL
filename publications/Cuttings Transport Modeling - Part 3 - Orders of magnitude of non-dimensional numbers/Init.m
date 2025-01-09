% Initialize

close all;
clear all;
clc;

% Add paths
p = pwd;
addpath([p(1) ':\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\alexabus']);
addpath([p(1) ':\Users\alexabus\OneDrive_NTNU\SWwd\MATLAB\Generic\MathWorks File Exchange']);
clear p;

% Color band
colorarray = 1:1:10; % Length depneds on length(PV)
N = length(colorarray); % Color definition: 
ColorSet = ColorBand(N); % from red to bpwd(1)lue
ColorSet = flip(ColorSet); % from blue to red
% Grey tones

% NTNU colors https://innsida.ntnu.no/wiki/-/wiki/Norsk/Farger+i+grafisk+profil
col_dp = [0 80 158]/255; % blue
col_an1 = [144 73 45]/255; % brown
col_an2 = [204, 189, 143]/255; % ligh brown
col_an3 = [223, 216, 197]/255; % lighter brown
lightblue = [121 162 206]/255; % light blue
purple = [173, 32, 142]/255;

% Unit conversion, multiply with the following conversion factors
in2m = 0.0254;
gpm2Lpm = 3.785411784;
gpm2m3ps = 0.00006309;
lbf100f22Pa = 0.4788026; % [lbf/100ft²] --> [Pa]
FannRPM2SR = 1.703; % Fann [rpm] --> [1/s]
FannDEG2lbf100f2 = 1.066;  % Fann [°] --> [lbf/100ft²]


% Ouriemi et al. (2010) flow pattern data
[~, ~, raw] = xlsread( [pwd '\Re_vs_Ar_Ouriemi_et_al_(2010).xlsx'],'Sheet1','A4:T51');
data = reshape([raw{:}],size(raw));
Ouriemi.FlatBed = [data(:,1), data(:,2)];
Ouriemi.FlatBed(~any(~isnan(Ouriemi.FlatBed), 2),:)=[];
Ouriemi.SmallDunes = [data(:,4), data(:,5)]; 
Ouriemi.SmallDunes(~any(~isnan(Ouriemi.SmallDunes), 2),:)=[];
Ouriemi.VortexDunes = [data(:,7), data(:,8)]; 
Ouriemi.VortexDunes(~any(~isnan(Ouriemi.VortexDunes), 2),:)=[];
Ouriemi.SinousDunes = [data(:,10), data(:,11)]; 
Ouriemi.SinousDunes(~any(~isnan(Ouriemi.SinousDunes), 2),:)=[];
Ouriemi.NoMotion = [data(:,13), data(:,14)]; 
Ouriemi.NoMotion(~any(~isnan(Ouriemi.NoMotion), 2),:)=[];
Ouriemi.Re_c = [data(:,16), data(:,17)]; 
Ouriemi.Re_c(~any(~isnan(Ouriemi.Re_c), 2),:)=[];
Ouriemi.theta_c = [data(:,19), data(:,20)]; 
Ouriemi.theta_c(~any(~isnan(Ouriemi.theta_c), 2),:)=[];
clearvars data raw;