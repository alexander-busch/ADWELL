% Creates folder structure required for TUI UDF compilation according to
% "5.3.1. Set Up the Directory Structure" in User Guide

clear all
clc

casepath = 'C:\Users\alexabus\IEPT1325_Local\Data\SW_working_directories\FLUENT\AdWell\Dune2D\';
currentFolder = pwd;
cd (casepath);

%% Create folders

mkdir libudf
mkdir libudf\src
mkdir libudf\win64
mkdir libudf\win64\2ddp




copyfile ('C:\Program Files\ANSYS Inc\v172\fluent\fluent17.2.0\src\udf\user_nt.udf', 'libudf\win64\2ddp');
copyfile ('C:\Program Files\ANSYS Inc\v172\fluent\fluent17.2.0\src\udf\makefile_nt.udf', 'libudf\win64\2ddp');

% Check why renaming to makefile does not work
oldname = [casepath 'libudf\win64\2ddp\makefile_nt.udf'];
newname = 'Makefile';
doscomand = strcat('rename', {' '}, oldname, {' '}, newname);


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