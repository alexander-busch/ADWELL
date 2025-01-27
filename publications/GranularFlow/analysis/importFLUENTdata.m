function [xcoordinate,ycoordinate,solidvof] = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [NODENUMBER,XCOORDINATE,YCOORDINATE,SOLIDVOF,VARNAME5] =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   [NODENUMBER,XCOORDINATE,YCOORDINATE,SOLIDVOF,VARNAME5] =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [nodenumber,xcoordinate,ycoordinate,solidvof,VarName5] = importfile('vof_init_cc.txt',2, 1068);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/11 16:43:07

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%10f%17f%17f%17f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Remove white space around all cell columns.
dataArray{5} = strtrim(dataArray{5});

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
nodenumber = dataArray{:, 1};
xcoordinate = dataArray{:, 2};
ycoordinate = dataArray{:, 3};
solidvof = dataArray{:, 4};
VarName5 = dataArray{:, 5};


