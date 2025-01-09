% If you plan to read the file with Microsoft® Notepad, use '\r\n' instead of '\n' to move to a new line.

clear all;
close all;
clc;

dim = '2D';
%dim = '3D';

load('Khatibi_et_al_2016.mat');

journalpath = ['C:\Users\alexabus\Data\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - ' dim];
resultspath = 'O:\\Data\\SWwd\\FLUENT\\AdWell\\Particle trajectories\\2017 - AB - 2D\\archive';
fileID = fopen([journalpath '\ParticleTrajectory_RunAllCases_' dim '.jou'],'w');

meshlist = {'case_refined'}; %, 'case_refined', 'case_with_traps', 'case_with_traps_refined'};

versionlist = {'_asis', '_new'};

i_max = size(Cases);


for cfdcase = 1:4
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, ['; ' num2str(cfdcase) ' \r\n']);
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    
    for version = 1:2 % '_asis' vs '_new'
        for i = 1:i_max(1)
       
        % Create name string 
        name = [ resultspath '\\' num2str(cfdcase) '\\cases_solved\\' Cases.fluid{i} '-cross_U' num2str(Cases.U(i)) '_d_p' num2str(Cases.d_p(i)) versionlist{version} ];    
        name_exp = [ resultspath '\\' num2str(cfdcase) '\\exports\\' Cases.fluid{i} '-cross_U' num2str(Cases.U(i)) '_d_p' num2str(Cases.d_p(i)) versionlist{version} ];    
        fprintf(fileID, ['/file/read-case-data "' name '.cas" OK\r\n']);
        fprintf(fileID, '/surface/line-surface line_xy_262mm 0.015 -0.005 0.262 0.020\r\n');
        fprintf(fileID, '/surface/line-surface line_xy_285mm 0.038 -0.005 0.285 0.020\r\n');
        fprintf(fileID, [ 'plot/plot y "' name_exp '_vx_xy_262mm.txt" n y 0 1 n n x-velocity line_xy_262mm ,\r\n']);
        fprintf(fileID, [ 'plot/plot y "' name_exp '_vx_xy_285mm.txt" n y 0 1 n n x-velocity line_xy_285mm ,\r\n']);
        fprintf(fileID, ['/file/write-case-data "' name '.cas" OK\r\n']);
        fprintf(fileID, '\r\n');

        end
    end
end

fclose(fileID);