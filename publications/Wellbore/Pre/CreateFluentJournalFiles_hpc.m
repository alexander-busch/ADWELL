% Create hpc array jobs journal files for multiphase (turbulent) annular flow
% Modify absolute paths of filepath_metadata & filepath_journal

clear all;
close all;
clc;

% Define case to simulate
%case2sim = '181026_1';

% Select 
% load('wellboredata.mat');
% case2sim = cases(strcmp(cases.ID,case2sim),:);

% Parameters
% if contains(case2sim.Turbulence,'SST')
%     turbmod = 'omega'; 
% else
%     turbmod = 'epsilon';
% end
turbmod = 'omega'; 
rpm = [0 100 200 300];
sweep = rpm;

% Path and filename of journalfile
filepath_journal = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Wellbore\hpc\journals_hpc';

% Purge
if exist(filepath_journal,'dir')
    rmdir(filepath_journal,'s');
end
mkdir(filepath_journal);

% Loop all dpdx
for ii = 1:length(sweep)

    fileID = fopen([filepath_journal '\run' num2str(ii) '.jou'],'w');

    % Write journal file
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; rpm sweep of multiphase annular flow\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; 26.10.2018, alexander.busch@alumni.ntnu.no\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; For hpc purposes: Reads generic case and sweeps pressure gradient\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Read case and change color scheme\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '/file/read-case case.cas\r\n');
    %fprintf(fileID, '/parallel/partition/print-active-partitions\r\n');
    fprintf(fileID, '/display/set/colors/color-scheme classic\r\n');
    fprintf(fileID, '/display/set/colors/background "white"\r\n');
    fprintf(fileID, '/display/set/colors/foreground "black"\r\n');
    fprintf(fileID, '/display/re-render\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Clean up working directory\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '!if exist libudf rmdir /s /q libudf\r\n');
    fprintf(fileID, '!if exist solution rmdir /s /q solution\r\n');
    fprintf(fileID, '!mkdir solution\r\n');
    fprintf(fileID, '!mkdir "solution/reports"\r\n');
    fprintf(fileID, '!mkdir "solution/exports"\r\n');
    fprintf(fileID, '!mkdir "solution/monitors"\r\n');
    fprintf(fileID, '!mkdir "solution/autosave"\r\n');
    fprintf(fileID, '!mkdir "solution/cfd-post"\r\n');
    fprintf(fileID, '!mkdir "solution/images"\r\n');
    
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Rotation\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, ['(rpsetvar ''rpm_pipe ' num2str(sweep(ii)) ' )\r\n']);
    fprintf(fileID, '(if (> (%%rpgetvar ''rpm_pipe) 0) (ti-menu-load-string (string-append "/define/boundary-conditions/set/wall wall_inner () mixture motion-bc y motion-bc-moving rotating y omega n " (number->string (/ (* 2 pi (%%rpgetvar ''rpm_pipe)) 60) ) " ai 1 aj 0 ak 0 q")) () )\r\n');
    
    
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Definition of time step (and mesh replacement)\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '(if (= (%%rpgetvar ''rpm_whirl) 0)\r\n');
    fprintf(fileID, '   (ti-menu-load-string (string-append "/solve/set/number-of-time-steps " (number->string ( round (/ (%%rpgetvar ''tmax) (%%rpgetvar ''dt) )))))\r\n');
    fprintf(fileID, '   (\r\n');
    fprintf(fileID, '   ; Definition of rounding\r\n');
    fprintf(fileID, '   (define (round-off z n) (let ((power (expt 10 n))) (/ (round (* power z)) power) ) )\r\n');
    fprintf(fileID, '   ; Define rounding precision, depends on time step\r\n');
    fprintf(fileID, '   (define precision (round (/ (log (/ 1 (%%rpgetvar ''dt))) (log 10)))) (display precision)\r\n');
    fprintf(fileID, '   ; Duration of one orbital cycle\r\n');
    fprintf(fileID, '   (define T (round-off (/ 1 (/ (abs (%%rpgetvar ''rpm_whirl)) 60)) precision)) (display T)\r\n');
    fprintf(fileID, '   ; Replace mesh after n orbital cycles\r\n');
    fprintf(fileID, '   (define n 4) (display n)\r\n');
    fprintf(fileID, '   ; Replace mesh after time_replacemesh seconds\r\n');
    fprintf(fileID, '   (define time_replacemesh (round-off (* (/ 1 (/ (abs (%%rpgetvar ''rpm_whirl)) 60)) n) precision)) (display time_replacemesh)\r\n');
    fprintf(fileID, '   ; Set number of timesteps\r\n');
    fprintf(fileID, '   (ti-menu-load-string (string-append "/solve/set/number-of-time-steps " (number->string ( round (/ time_replacemesh (%%rpgetvar ''dt) )))))\r\n');
    fprintf(fileID, '   ; Number of required inner loops to achieve total flow time\r\n');
    fprintf(fileID, '   (define jMax (round (/ (%%rpgetvar ''tmax) time_replacemesh))) (display jmax)\r\n');
    fprintf(fileID, '   )\r\n');
    fprintf(fileID, ')\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
%     fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
%     fprintf(fileID, '; Reset dpdx or mdot\r\n');
%     fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
%     fprintf(fileID, '\r\n');
%     if strcmp(turbmod,'epsilon')
%         fprintf(fileID, ['/define/periodic-conditions/massflow-rate-specification ' num2str(m_dot(ii)) ' ' num2str(dpdx(ii)) ' 0.5 5 1 0 0\r\n']);
%         %fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 3 y y "solution/monitors/mdot=' num2str(m_dot(ii)) '_mass-flow-rate.txt" 1\r\n']);
%         %fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 4 y y "solution/monitors/mdot=' num2str(m_dot(ii)) '_bulk-velocity.txt" 1\r\n']);
%         fprintf(fileID, '/solve/monitors/statistic/monitors y\r\n'); 
%         fprintf(fileID, '/solve/monitors/statistic/print y\r\n'); 
%         fprintf(fileID, '/solve/monitors/statistic/plot y\r\n'); 
%         fprintf(fileID, '/solve/monitors/statistic/write y\r\n'); 
%         fprintf(fileID, '/solve/monitors/statistic/window 5\r\n');
%     else
%         fprintf(fileID, ['/define/periodic-conditions/pressure-gradient-specification ' num2str(dpdx(ii)) ' 1 0 0\r\n']);
%         %fprintf(fileID, ['/solve/monitors/surface/set-monitor mass-flow-rate "Mass Flow Rate" plane_yz () y 2 y y "solution/monitors/dpdx=' num2str(dpdx(ii)) '_mass-flow-rate.txt" 1\r\n']);
%         %fprintf(fileID, ['/solve/monitors/surface/set-monitor bulk-velocity "Area-Weighted Average" x-velocity plane_yz () y 3 y y "solution/monitors/dpdx=' num2str(dpdx(ii)) '_bulk-velocity.txt" 1\r\n']);
%     end
%     fprintf(fileID, '\r\n');
%     fprintf(fileID, '\r\n');
    
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Initialize\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, ['/solve/initialize/set-defaults/fluid/x-velocity 0\r\n']); %' num2str(vel(ii)) '
    fprintf(fileID, '/solve/initialize/set-defaults/fluid/k 1\r\n'); % num2str(k(ii)) '\r\n']);
    fprintf(fileID, ['/solve/initialize/set-defaults/fluid/' turbmod ' 1 \r\n']); % ' num2str(epsilon(ii)) '\r\n']);
    fprintf(fileID, ['/solve/initialize/set-defaults/solid/x-velocity 0\r\n']); %' num2str(vel(ii)) '
    fprintf(fileID, '/solve/initialize/initialize-flow ok\r\n');
    fprintf(fileID, ['/file/write-case-data "case_init.cas" OK\r\n']);
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '; Solve\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; Define pressure gradient\r\n');
    fprintf(fileID, '(define dpdx 500) (display dpdx)\r\n');
    fprintf(fileID, '(define dpdx_delta 500) (display dpdx)\r\n');
    fprintf(fileID, '(define dpdx_max (+ dpdx (* dpdx_delta 6))) (display dpdx_max)\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; Loop all pressure gradients\r\n');
    fprintf(fileID, '(do ((i dpdx (+ i dpdx_delta)))\r\n');
    fprintf(fileID, '	((> i dpdx_max))\r\n');
    fprintf(fileID, '       (display dpdx)\r\n');
    fprintf(fileID, '       ; Set RP variabble for pressure gradient\r\n');
    fprintf(fileID, '       (rpsetvar ''mixture/dpdx (* dpdx -1))\r\n');
    fprintf(fileID, '       ; Update periodic BC for pressure gradient\r\n');
    fprintf(fileID, '       (ti-menu-load-string (string-append "/define/periodic-conditions/pressure-gradient-specification " (number->string (%%rpgetvar ''mixture/dpdx)) " 1 0 0"))\r\n');
    fprintf(fileID, '       (if (= (%%rpgetvar ''rpm_whirl) 0)\r\n');
    fprintf(fileID, '           (\r\n');
    fprintf(fileID, '           ; Solve\r\n');
    fprintf(fileID, '           (ti-menu-load-string "/solve/dual-time-iterate , ,")\r\n');  
    fprintf(fileID, '           )\r\n');
	fprintf(fileID, '           (\r\n');
	fprintf(fileID, '           ; Loop all mesh replacement intervalls and solve per interval, mesh is replaced after n intervals to avoid mesh degration\r\n');
	fprintf(fileID, '           (ti-menu-load-string "/solve/dual-time-iterate , ,")\r\n');  
    fprintf(fileID, '           (do ((j 1 (+ j 1))) ((= j (- jMax 1))) (ti-menu-load-string "/mesh/replace case.cas n") (ti-menu-load-string "/solve/dual-time-iterate , , y y y y y y y") )\r\n');
    fprintf(fileID, '           )\r\n');
    fprintf(fileID, '       )\r\n');
    fprintf(fileID, '       ; Write solved cases to disk\r\n');
    fprintf(fileID, '       (ti-menu-load-string (string-append "/file/write-case-data " (string-append "solution/dpdx" (number->string (%%rpgetvar ''mixture/dpdx)) "_solved" ) " OK"))\r\n');
    fprintf(fileID, '       ; Increase pressure gradient\r\n');
    fprintf(fileID, '       (set! dpdx (+ dpdx dpdx_delta))\r\n');
    fprintf(fileID, '       )\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '/report/summary y "solution/reports/casesetup.txt"\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '\r\n');
    
    fclose(fileID);
end

