% Compare various drag laws

%% Initialize

% Clean-up
clear all;
close all;
clc;

% Settings
format long;

% Folders
addpath 'E:\OneDrive - NTNU\SWwd\MATLAB\Generic\m-files';
addpath([pwd '\DragLaws']);
addpath([pwd '\Utilities'])

% Parameters
Definitions;
Parameters;


% Create cases table

% Create
Cases = table;

% Starting row index
row = 1;

% Loop all fluids, particle diameters and build case table
for i = 1:length(fluidlist)
    for j = 1:length(d_p)
        
        % Build Cases table
        Cases.index(row,1) = row;
        Cases.fluid(row,1) = string(fluidlist{i});
        Cases.d_p(row,1) = d_p(j);
        
        % Increase row index
        row = row + 1;
        
    end % of j
end % of i


% Size opf cases
i_max = size(Cases);


%% Compute drag

% Loop all drag laws
for i = 1:length(draglawlist)
    
    % Get current drag law
    draglaw = draglawlist{i};
    
    % Create table for solution of current drag law
    eval([draglawlist{i,2} ' = table']);
    
    % Loop all cases
    for j = 1:i_max(1)
        
        % Determine row index of Cross table for current fluid
        k = find(Cross.Fluid==Cases.fluid(j));
    
        IterateSettlingVelocity;
        
    end % of j
    
end % of i


%% Compare drag

% Sort SN1935 in ascending order by Re_p
[SN1935, idx] = sortrows(SN1935,{'Re_p'},{'ascend'});

% Sort all othe tables in order by idx
A1976 = A1976(idx,:);
CU1980 = CU1980(idx,:);
C1999 = C1999(idx,:);
R2004 = R2004(idx,:);
S2007 = S2007(idx,:);

CreateFigure('Drag law comparison', 'Re_p [-]', 'c_D_-_i/d_D_-_S_N_1_9_3_5 [-]');

    Legend = cell(length(draglawlist),1);
    
    plot(SN1935.Re_p2,SN1935.c_D./SN1935.c_D,'k-o');
    plot(A1976.Re_p2,A1976.c_D./SN1935.c_D,'m-o');
    plot(CU1980.Re_p2,CU1980.c_D./SN1935.c_D,'b-o');
    plot(C1999.Re_p2,C1999.c_D./SN1935.c_D,'y-o');
    plot(R2004.Re_p2,R2004.c_D./SN1935.c_D,'r-o');
    plot(S2007.Re_p2,S2007.c_D./SN1935.c_D,'g-o');
    
    for i=1:length(draglawlist)
        Legend{i} = draglawlist{i};
    end 
    
    set(gca,...
        'XScale','log',...
        'YScale','lin',...
        'YLim',[0.01 2]);
    
    legend(Legend, 'Location','northeast');

CreateFigure('Drag law comparison', 'Re_p [-]', 'c_D_i [-]');
    
    Legend = cell(length(draglawlist),1);
    
    plot(SN1935.Re_p2,SN1935.c_D,'k-o');
    plot(A1976.Re_p2,A1976.c_D,'m-o');
    plot(CU1980.Re_p2,CU1980.c_D,'b-o');
    plot(C1999.Re_p2,C1999.c_D,'y-o');
    plot(R2004.Re_p2,R2004.c_D,'r-o');
    plot(S2007.Re_p2,S2007.c_D,'g-o');
    
    
    for i=1:length(draglawlist)
        Legend{i} = draglawlist{i};
    end
    
    legend(Legend, 'Location','northeast');
    
    set(gca,...
        'XScale','log',...
        'YScale','log');
    
    
    
%% Compare velocity

CreateFigure('Settling velocity comparison', 'Re_p [-]', 'v_s_e_t_-_i / v_s_e_t_-_S_N_1_9_3_5 [m/s]');

    Legend = cell(length(draglawlist),1);
   
    plot(SN1935.Re_p2,SN1935.v_set./SN1935.v_set,'k-o');
    plot(A1976.Re_p2,A1976.v_set./SN1935.v_set,'m-o');
    plot(CU1980.Re_p2,CU1980.v_set./SN1935.v_set,'b-o');
    plot(C1999.Re_p2,C1999.v_set./SN1935.v_set,'y-o');
    plot(R2004.Re_p2,R2004.v_set./SN1935.v_set,'r-o');
    plot(S2007.Re_p2,S2007.v_set./SN1935.v_set,'g-o');
        
    for i=1:length(draglawlist)
        Legend{i} = draglawlist{i};
    end 
    
    legend(Legend, 'Location','northeast');
    
    set(gca,...
        'XScale','log',...
        'YScale','lin',...
        'YLim',[0.5 2]);
