function y = importheaders(filename)
% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\alexanderb\Documents\MATLAB\PAC4.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2015/12/08 12:20:58

%% Initialize variables.
delimiter = ',';
startRow = 30;
endRow = 31;

%% Format string for each line of text:
%   column2: text (%s)
%	column3: text (%s)
%   column4: text (%s)
%	column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%s%s%s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
y = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;