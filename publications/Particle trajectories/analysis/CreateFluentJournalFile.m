% If you plan to read the file with Microsoft® Notepad, use '\r\n' instead of '\n' to move to a new line.

clear all;
close all;
clc;


%% Definitions


% CFD problems
problemlist = {'B'}; 
% A = Long with traps
% B = Short with traps
% C = Short without traps

% Slip velocity at wall
slip = {'noslip'};
    
% Spatial dimensions
if strcmp(problemlist,'A')
    dim = '3D';
else
    dim = '2D';
end

% Rotation
rotation = 'norotation'; % rotation

% Model version: Fluent as is vs. new shear rate/drag law model
versionlist = {'_asis', '_new'};

% Shear rate particle model
model_SR_p = 1;
% 1 = v_r / d_p; e.g. Acharya (1976)
% 2 = 2 v_r / d_p; e.g. Chhabra & Uhlherr (1980) 
% 3 = sqrt(6) v_r / d_p; Renaud et al. (2004)
% 4 = tbd 

% Shear rate model
model_SR = 1;
% 1 = Vectorsum
% 2 = based on unitvectors +
% 3 = based on unitvectors -
% 4 = tbd 

% Rheology model
model_R = '"Cross"';
% "Cross"
% "Carreau"


% Drag law model
model_c_D = '"SN"';
% "SN" = Schiller-Naumann (1935)
% "A" = Acharya (1976) - requires SR =1
% "CU" = Chhabra & Uhlherr (1980)-  requires SR =1
% "FS" = Fluents spherical drag law


% Path and filename of journalfile
journalpath = ['O:\OneDrive_NTNU\SWwd\FLUENT\AdWell\Particle trajectories\2017 - AB - ' dim];
fileID = fopen([journalpath '\ParticleTrajectory_RunAllCases_' dim '.jou'],'w');


% Experimental data
load('Khatibi_et_al_2016.mat');

% Rheometric data
load Cross;
load Carreau;

% Length of caseslist
i_max = size(Cases);


%% Write journal file

% Loop all problems
for problem = 1:length(problemlist)
    fprintf(fileID, '\r\n');
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, ['; ' problemlist{problem} ' - Working directory clean up & load case\r\n']);
    fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, ['!if exist archive\\case_' dim '_' problemlist{problem} ' rmdir /s /q archive\\case_' dim '_' problemlist{problem} '\r\n']);
    fprintf(fileID, '\r\n');
    
    % Clean up
    fprintf(fileID, '!if exist cases rmdir /s /q cases\r\n');
    fprintf(fileID, '!mkdir cases\r\n');
    fprintf(fileID, '!if exist cases_solved rmdir /s /q cases_solved\r\n');
    fprintf(fileID, '!mkdir cases_solved\r\n');
    fprintf(fileID, '!if exist exports rmdir /s /q exports\r\n');
    fprintf(fileID, '!mkdir exports\r\n');
    fprintf(fileID, '!if exist images rmdir /s /q images\r\n');
    fprintf(fileID, '!mkdir images\r\n');
    fprintf(fileID, '!if exist DPM_x.txt del DPM_x.txt\r\n');
    fprintf(fileID, '!if exist DPM_y.txt del DPM_y.txt\r\n');
    fprintf(fileID, '!if exist DPM_rep.txt del DPM_rep.txt\r\n'); 
    fprintf(fileID, '\r\n');
    
    % Read & compile
	if strcmp(dim,'2D')
        fprintf(fileID, ['/file/read-case-data "case_' dim '_' problemlist{problem} '.cas" OK\r\n']);
    else
        fprintf(fileID, ['/file/read-case-data "case_' dim '_A.cas" OK\r\n']);
    end
                        
    fprintf(fileID, '/define/user-defined/compiled-functions load "libudf"\r\n');
    fprintf(fileID, '\r\n');
    
    % Reset convergence criteria depending on dimension
    if strcmp(dim,'2D')
        fprintf(fileID, '/solve/monitors/residual/convergence-criteria 1e-5 1e-5 1e-5\r\n');
    else
        fprintf(fileID, '/solve/monitors/residual/convergence-criteria 1e-5 1e-5 1e-5 1e-5\r\n');
    end
    
    fprintf(fileID, '\r\n');
    
    % Set shear rate particle model
    fprintf(fileID, '(rpsetvar ''model/sr_p %u)\r\n', model_SR_p);
        
    % Set shear rate model
    fprintf(fileID, '(rpsetvar ''model/sr %u)\r\n', model_SR);

    % Rheology model
    fprintf(fileID, '(rpsetvar ''model/rheology %s)\r\n', model_R);

    % Drag law model
    fprintf(fileID, '(rpsetvar ''model/drag %s)\r\n', model_c_D);
    
    fprintf(fileID, '\r\n');
        
% Loop all model versions
for version = 1:2 % '_asis' vs '_new'
    
    % Loop all CFD cases
    for i = 1:i_max(1)
        
        % Full name of current case
        currentcase = ['Problem=' problemlist{problem} '_Version=' versionlist{version} '_' Cases.fluid{i,1} '_d_p=' num2str(Cases.d_p(i,1)) '_U=' num2str(Cases.U(i,1))];
        % fprintf(fileID, ['(rpsetvar ''currentcase ' currentcase ')\r\n']);
        
        fprintf(fileID, '\r\n');
        fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
        fprintf(fileID, ['; ' currentcase '\r\n']);
        fprintf(fileID, '; --------------------------------------------------------------------------------------------\r\n');
        fprintf(fileID, '\r\n');

        % Drag law depending on version
        fprintf(fileID, '/define/models/dpm/injections/properties/set/pick-injections-to-set n injection , n\r\n');
        if version == 1
            % fprintf(fileID, '/define/models/dpm/injections/properties/set/physical-models/drag-parameters spherical\r\n');
            fprintf(fileID, '/define/models/dpm/injections/properties/set/physical-models/drag-parameters drag_force_n::libudf\r\n');
        else
            fprintf(fileID, '/define/models/dpm/injections/properties/set/physical-models/drag-parameters drag_force_st::libudf\r\n');
        end
        fprintf(fileID, '\r\n');
        
        % Fluid
        
        % Laminar flow model
        fprintf(fileID, '/define/models/viscous/laminar y\r\n');
        % fprintf(fileID, '/define/models/viscous/kw-sst y\r\n');
        % fprintf(fileID, '/define/models/viscous/kw-low-re-correction y\r\n');
        % fprintf(fileID, '/solve/set/discretization-scheme k 1\r\n');
        % fprintf(fileID, '/solve/set/discretization-scheme omega 1\r\n');
        % fprintf(fileID, '/solve/monitors/residual/convergence-criteria 1e-7 1e-7 1e-7 1e-7 1e-7\r\n');
        
        % H2O
        if strcmp(Cases.fluid{i,1},'H2O')
            
            % Fluid variable for UDF
            fprintf(fileID, '(rpsetvar ''fluid 1)\r\n');
            
            fprintf(fileID, '(rpsetvar ''cross/lambda 0.0)\r\n');
            fprintf(fileID, '(rpsetvar ''cross/n 1.0)\r\n');
            fprintf(fileID, '(rpsetvar ''cross/mu_0 0.001002)\r\n');
            fprintf(fileID, '(rpsetvar ''cross/mu_inf 0.001002)\r\n');
            
            fprintf(fileID, '(rpsetvar ''carreau/lambda 0.0)\r\n');
            fprintf(fileID, '(rpsetvar ''carreau/n 1.0)\r\n');
            fprintf(fileID, '(rpsetvar ''carreau/mu_0 0.001002)\r\n');
            fprintf(fileID, '(rpsetvar ''carreau/mu_inf 0.001002)\r\n');
            
            % BC
            if strcmp(dim,'2D')
                fprintf(fileID, '/define/boundary-conditions fluid , y h2o , , , , , , , , , ,\r\n');
            else
                fprintf(fileID, '/define/boundary-conditions fluid , y h2o , , , , , , , , , , , , , , , , , , ,\r\n');
            end
             
        % PAC
        else
            if strcmp(Cases.fluid{i,1},'PAC2')
            
                % Fluid variable for UDF
                fprintf(fileID, '(rpsetvar ''fluid 2)\r\n');

                % Cross model parameters (from own curve fit of Khatibi et al.
                % (2016) data
                fprintf(fileID, '(rpsetvar ''cross/lambda %2.5f)\r\n', Cross.lambda{2});
                fprintf(fileID, '(rpsetvar ''cross/n %2.5f)\r\n', Cross.n{2});
                fprintf(fileID, '(rpsetvar ''cross/mu_0 %2.5f)\r\n', Cross.mu_0{2});
                fprintf(fileID, '(rpsetvar ''cross/mu_inf %2.5f)\r\n', Cross.mu_inf{2});

                % Carreau model parameters (as given in Khatibi et al. (2016)
                % To be corrected with own curve fit as muinf is incorrect
                fprintf(fileID, '(rpsetvar ''carreau/lambda %2.5f)\r\n', Carreau.lambda{2});
                fprintf(fileID, '(rpsetvar ''carreau/n %2.5f)\r\n', Carreau.n{2});
                fprintf(fileID, '(rpsetvar ''carreau/mu_0 %2.5f)\r\n', Carreau.mu_0{2});
                fprintf(fileID, '(rpsetvar ''carreau/mu_inf %2.5f)\r\n', Carreau.mu_inf{2});
                        
            elseif strcmp(Cases.fluid{i,1},'PAC4')
                        
                % Fluid variable for UDF
                fprintf(fileID, '(rpsetvar ''fluid 3)\r\n');

                % Cross model parameters (from own curve fit of Khatibi et al.
                % (2016) data
                fprintf(fileID, '(rpsetvar ''cross/lambda %2.5f)\r\n', Cross.lambda{3});
                fprintf(fileID, '(rpsetvar ''cross/n %2.5f)\r\n', Cross.n{3});
                fprintf(fileID, '(rpsetvar ''cross/mu_0 %2.5f)\r\n', Cross.mu_0{3});
                fprintf(fileID, '(rpsetvar ''cross/mu_inf %2.5f)\r\n', Cross.mu_inf{3});

                % Carreau model parameters (as given in Khatibi et al. (2016)
                % To be corrected with own curve fit as muinf is incorrect
                fprintf(fileID, '(rpsetvar ''carreau/lambda %2.5f)\r\n', Carreau.lambda{3});
                fprintf(fileID, '(rpsetvar ''carreau/n %2.5f)\r\n', Carreau.n{3});
                fprintf(fileID, '(rpsetvar ''carreau/mu_0 %2.5f)\r\n', Carreau.mu_0{3});
                fprintf(fileID, '(rpsetvar ''carreau/mu_inf %2.5f)\r\n', Carreau.mu_inf{3});
                
            else
            end
            
            % Specification of viscosity model to be used
            if strcmp(model_R, '"Cross"')
                
                % Fluents Cross model (mu_inf = 0)
                % fprintf(fileID, '/define/materials/change-create pac-cross pac-cross n n n y cross shear-rate-dependent (%%rpgetvar ''cross/mu_0) (%%rpgetvar ''cross/n) (%%rpgetvar ''cross/lambda) n n n n\r\n');
                
                % Fluents UDF Cross model (mu_inf = 0.00102)
                fprintf(fileID, '/define/materials/change-create pac-cross pac-cross n n n y user-defined "rheology_cross::libudf" n n n\r\n');
                
            elseif strcmp(model_R, '"Carreau"')
                
                % Fluents Carreau model
                fprintf(fileID, '/define/materials/change-create pac-carreau pac-carreau n n n y carreau shear-rate-dependent (%%rpgetvar ''carreau/lambda) (%%rpgetvar ''carreau/n) (%%rpgetvar ''carreau/mu_0) (%%rpgetvar ''carreau/mu_inf) n n n\r\n');
                
            else
            end
            
            % BC
            if strcmp(dim,'2D')
                if strcmp(model_R, '"Cross"')
                    
                    fprintf(fileID, '/define/boundary-conditions fluid , y pac-cross , , , , , , , , , , ,\r\n');
                    
                elseif strcmp(model_R, '"Carreau"')
                    
                    fprintf(fileID, '/define/boundary-conditions fluid , y pac-carreau , , , , , , , , , ,\r\n');
                    
                else
                end

            else
                if strcmp(model_R, '"Cross"')
                    
                    fprintf(fileID, '/define/boundary-conditions fluid , y pac-cross , , , , , , , , , , , , , , , , , , ,\r\n');
                    
                elseif strcmp(model_R, '"Carreau"')
                    
                    fprintf(fileID, '/define/boundary-conditions fluid , y pac-carreau , , , , , , , , , , , , , , , , , ,\r\n');
                    
                else
                end
            end
            
        end
        
        fprintf(fileID, '\r\n');

        % Particle diameter
        fprintf(fileID, '(rpsetvar ''d_p %2.5f)\r\n', Cases.d_p(i,1));
        fprintf(fileID, '\r\n');
        
        % Set particle injection properties. Actually this is a function of
        % problemlist, B, C = 0.008 and A = 0.255
        % USe 0.019999999 instead of 0.02 in order to ensure that particles
        % are injected in the domain.
        if strcmp(dim,'2D')
            if strcmp(rotation,'rotation')
                fprintf(fileID, '/define/models/dpm/injections/set-injection-properties injection , n n n n n n n 0.008 0.019999999 %2.7f %2.7f 0 (%%rpgetvar ''d_p) 0.0094\r\n', EXP.vx0(i,1), EXP.vy0(i,1));
            else
                fprintf(fileID, '/define/models/dpm/injections/set-injection-properties injection , n n n n n n n 0.008 0.019999999 %2.7f %2.7f (%%rpgetvar ''d_p) 0.0094\r\n', EXP.vx0(i,1), EXP.vy0(i,1));
            end
        else
            if strcmp(rotation,'rotation')
                fprintf(fileID, '/define/models/dpm/injections/set-injection-properties injection , n n n n n n n 0.255 0.019999999 0.000 %2.7f %2.7f 0 0 0 0 (%%rpgetvar ''d_p) 0.0094\r\n', EXP.vx0(i,1), EXP.vy0(i,1));
            else
                fprintf(fileID, '/define/models/dpm/injections/set-injection-properties injection , n n n n n n n 0.255 0.019999999 0.000 %2.7f %2.7f 0 (%%rpgetvar ''d_p) 0.0094\r\n', EXP.vx0(i,1), EXP.vy0(i,1));
            end
        end       
        fprintf(fileID, '\r\n');

        % Fluid velocity
        fprintf(fileID, '(rpsetvar ''u_bulk %2.4f)\r\n', Cases.U(i,1));
        if strcmp(dim,'2D')
            % fprintf(fileID, '/define/boundary-conditions/set/velocity-inlet inlet , vmag y y "udf" "exp_vel_pro::libudf" q\r\n');
            fprintf(fileID, '/define/boundary-conditions/set/velocity-inlet inlet , vmag y y "udf" "par_vel_pro::libudf" q\r\n');
        else
            fprintf(fileID, '/define/boundary-conditions/set/velocity-inlet inlet , vmag n (%%rpgetvar ''u_bulk) q\r\n'); % Bulk velocity
        end

        
        % Get slip velocities at top wall        
        % Slip velocities at bottom wall are x-times slip velocity at top wall
        if strcmp(Cases.fluid{i,1},'H2O') % H2O - No slip
            fprintf(fileID, '(rpsetvar ''u_slip_top %2.5f)\r\n', 2*EXP.u_slip(i,1));
            fprintf(fileID, '(rpsetvar ''u_slip_bottom (* 1.1 (%%rpgetvar ''u_slip_top)))\r\n');
        elseif strcmp(Cases.fluid{i,1},'PAC2') % PAC - Slip
            fprintf(fileID, '(rpsetvar ''u_slip_top %2.5f)\r\n', EXP.u_slip(i,1));
            fprintf(fileID, '(rpsetvar ''u_slip_bottom (* 5.5 (%%rpgetvar ''u_slip_top)))\r\n');
        elseif strcmp(Cases.fluid{i,1},'PAC4') % PAC - Slip
            fprintf(fileID, '(rpsetvar ''u_slip_top %2.5f)\r\n', 0.4*EXP.u_slip(i,1));
            fprintf(fileID, '(rpsetvar ''u_slip_bottom (* 5.5 (%%rpgetvar ''u_slip_top)))\r\n');
        else
        end
        fprintf(fileID, '\r\n');
        
        % Assign slip velocities
        if strcmp(slip,'slip')
            fprintf(fileID, '/define/boundary-conditions/set/wall wall_top () vmag n (%%rpgetvar ''u_slip_top) q\r\n');
            fprintf(fileID, '\r\n');
            if strcmp(problemlist,'C')
                fprintf(fileID, '/define/boundary-conditions/set/wall wall_bottom () vmag n (%%rpgetvar ''u_slip_bottom) q\r\n');
            else
                fprintf(fileID, '/define/boundary-conditions/set/wall wall_bottom_slip () vmag n (%%rpgetvar ''u_slip_bottom) q\r\n');
            end
        else
            fprintf(fileID, '/define/boundary-conditions/set/wall wall_top () vmag n 0.0 q\r\n');
            fprintf(fileID, '\r\n');
            if strcmp(problemlist,'C')
                fprintf(fileID, '/define/boundary-conditions/set/wall wall_bottom () vmag n 0.0 q\r\n');
            else
                fprintf(fileID, '/define/boundary-conditions/set/wall wall_bottom_slip () vmag n 0.0 q\r\n');
            end
        end
        fprintf(fileID, '\r\n');

        % Create name string
        name = [ Cases.fluid{i} '-cross_U' num2str(Cases.U(i)) '_d_p' num2str(Cases.d_p(i)) versionlist{version} ];
        
        % Initialize & solve
        fprintf(fileID, '/solve/initialize/initialize-flow y\r\n');
        fprintf(fileID, [ '/file/write-case-data "cases/' name '.cas" OK\r\n']);
        fprintf(fileID, '\r\n');
        if strcmp(dim,'2D')
        	fprintf(fileID, '/solve/iterate 600\r\n');
        else
        	fprintf(fileID, '/solve/iterate 1000\r\n');
        end
        fprintf(fileID, '\r\n');
                        
        fprintf(fileID, [ '/file/write-case-data "cases_solved/' name '.cas" OK\r\n']);
        fprintf(fileID, '\r\n');

        % Post
        fprintf(fileID, '/display/set-window 2\r\n');
        if strcmp(dim,'2D')
        	
        else
            fprintf(fileID, '/display/set/contours/surfaces z=0 ()\r\n');
        end
       
        fprintf(fileID, '/views/restore-view xy\r\n');
        fprintf(fileID, '/display/contour/velocity-magnitude , ,\r\n');
        fprintf(fileID, '/display/set/overlays y\r\n');
        fprintf(fileID, '/display/particle-tracks particle-tracks particle-reynolds-number , , , ,\r\n');
        fprintf(fileID, '/display/particle-tracks particle-tracks particle-reynolds-number , , , ,\r\n');
        fprintf(fileID, '\r\n');

        fprintf(fileID, [ '/display/save-picture "images/' name '.png" OK\r\n']);
        fprintf(fileID, '\r\n');
                
        fprintf(fileID, '/display/set-window 3\r\n');
        
        if strcmp(dim,'2D')
        	fprintf(fileID, [ 'plot/plot y "exports/' name '_vx_xy.txt" n y 0 1 n n x-velocity line_xy ,\r\n']);
        else
            fprintf(fileID, [ 'plot/plot y "exports/' name '_vx_xy.txt" n y 0 1 0 n n x-velocity line_xy ,\r\n']);
            fprintf(fileID, [ 'plot/plot y "exports/' name '_vx_xz.txt" n y 0 0 1 n n x-velocity line_xz ,\r\n']);
        end
        fprintf(fileID, '\r\n');
        fprintf(fileID, '\r\n');
        
        % Particle trajectory exports
        fprintf(fileID, '; Particle trajectory exports\r\n');
        fprintf(fileID, '; These are conducted via GUI as this yields better data quality, see xlsx file in /trials\r\n');
                
        if strcmp(dim,'2D')
        	fprintf(fileID, '(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Results|Graphics|Particle Tracks"))\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Results|Graphics|Particle Tracks"))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "NavigationPane*List_Tree1")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)" #f)\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)" #t)\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)" #t)\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 3))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList4(X Axis Function)" ''( 0))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList4(X Axis Function)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*List5(Release from Injections)" ''())\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*List5(Release from Injections)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*List5(Release from Injections)" ''( 0))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*List5(Release from Injections)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_x.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 4))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_y.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 22))\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_rep.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, '(cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton2(Cancel)")\r\n');
        else
            fprintf(fileID, ' (cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Results|Graphics|Particle Tracks"))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Results|Graphics|Particle Tracks"))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "NavigationPane*List_Tree1")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)" #f)\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)" #f)\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)" #t)\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton4(XY Plot)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-toggle-button2 "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)" #t)\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Table1*Frame1*Frame1*ToggleBox1(Options)*CheckButton5(Write to File)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 3))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList4(X Axis Function)" ''( 0))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList4(X Axis Function)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*List5(Release from Injections)" ''())\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*List5(Release from Injections)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*List5(Release from Injections)" ''( 0))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*List5(Release from Injections)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_x.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 4))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_y.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-list-selections "Particle Tracks*Frame2*Table1*DropDownList2" ''( 22))\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*Frame2*Table1*DropDownList2")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton1(OK)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-set-file-dialog-entries "Select File" ''( "DPM_rep.txt") "Particle Tracks Files (*.xy)")\r\n');
            fprintf(fileID, ' (cx-gui-do cx-activate-item "Particle Tracks*PanelButtons*PushButton2(Cancel)")\r\n');
        end
        fprintf(fileID, '\r\n');
        
        fprintf(fileID, ['!move /y DPM_x.txt exports\\' name '_x.txt\r\n']);
        fprintf(fileID, ['!move /y DPM_y.txt exports\\' name '_y.txt\r\n']);
        fprintf(fileID, ['!move /y DPM_rep.txt exports\\' name '_rep.txt\r\n']);
         
%         fprintf(fileID, [ '/display/particle-tracks/plot-write-xy-plot x-coordinate , , time , , y "exports/' name '_x.txt" OK\r\n']);
%         fprintf(fileID, [ '/display/particle-tracks/plot-write-xy-plot y-coordinate , , time , , y "exports/' name '_y.txt" OK\r\n']);
%         fprintf(fileID, [ '/display/particle-tracks/plot-write-xy-plot particle-reynolds-number , , time , , y "exports/' name '_rep.txt" OK\r\n']);
        fprintf(fileID, '\r\n');
        
        fprintf(fileID, [ '/file/write-case-data "cases_solved/' name '.cas" OK\r\n']);
        fprintf(fileID, '\r\n');
    end
end

fprintf(fileID, ['!mkdir archive\\case_' dim '_' problemlist{problem} '\r\n']);
fprintf(fileID, ['!move /y cases archive\\case_' dim '_' problemlist{problem} '\\cases\r\n']);
fprintf(fileID, ['!move /y cases_solved archive\\case_' dim '_' problemlist{problem} '\\cases_solved\r\n']);
fprintf(fileID, ['!move /y exports archive\\case_' dim '_' problemlist{problem} '\\exports\r\n']);
fprintf(fileID, ['!move /y images archive\\case_' dim '_' problemlist{problem} '\\images\r\n']);

fprintf(fileID, ['!robocopy /e libudf archive\\case_' dim '_' problemlist{problem} '\\cases\\libudf\r\n']);
fprintf(fileID, ['!robocopy /e libudf archive\\case_' dim '_' problemlist{problem} '\\cases_solved\\libudf\r\n']);
    
end

fclose(fileID);