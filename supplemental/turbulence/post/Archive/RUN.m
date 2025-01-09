close all;
clear all;
clc;

read = 'm_dot'; % 'm_dot'

geo = 'p'; % 'a' 'p'
% fluid = 'h2o';
%fluid = 'pl';
fluid = 'cross';
dpdx = [0.001; 0.0025; 0.005; 0.01; 0.013; 0.01447; 0.018; 0.02; 0.03; 0.04];




filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\archive';

cases = readtable([filepath '\Cases.xlsx']);

% Paths
currentpath = pwd;
parentpath = cd(cd('..'));
archivepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Sand Pile\17.2\archive';

cd(filepath)
cd(currentpath)
% Purge and create folder for CFD post exports
if exist([parentpath '\cfdpostexports'],'dir')


% Loop all table rows
for aa=1:height(cases)
    % Debugging aa = 5
   
    for bb=1:length(dpdx)
        
        
        
        
        if strcmp(cases.Specification{aa}, 'dpdx')
        	filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_dpdl=' num2str(dpdx(bb)) '_mass-flow-rate.txt'];
        else
            % Estimate mdot based on dpdx (Blasius)
            U=(dpdx(bb)^4*cases.d_h(aa)^5/(1/16*0.316^4*cases.scaledK_PL(aa)*cases.scaledrho(aa)^3))^(1/7);
            A=pi/4*(cases.d_o(aa)^2-cases.d_i(aa)^2);
            Q=U*A;
            mdot(bb)=Q/cases.scaledrho(aa);
    
            filename = ['\' cases.Geometry{aa}(1) '_' cases.Fluid{aa} '-scaled_mdot=' num2str(m_dot(bb)) '_dpdl.txt'];
        end
                
    end
        
    
    % Change directory 
    cd([filepath '\' cases.ID{aa}]);
    % Check folder structure and assemble correct path to results
    if exist('monitors')
        fullpath = [filepath '\' cases.ID{aa}'\monitors'  filename];
    else
        cd
    
    end
    
    
    
    

    
    
    
    
    
    
    
end



name = [filepath '\' cases.ID{aa}]

name = 'test'
7==exist(name,'dir')


raw{:,2}


% [~, ~, raw] = xlsread(filepath_metadata,'CFD','BD93:BJ102');
data = reshape([raw{:}],size(raw));
dpdl = data(:,1);
m_dot = data(:,7);

for ii=1:length(dpdl)

    if strcmp(read, 'm_dot')
        filename = ['\' geo '_' fluid '-scaled_dpdl=' num2str(dpdl(ii)) '_mass-flow-rate.txt'];
%         fullpath = [filepath filename];
        fullpath = [filepath '\' num2str(ii) '\solution\monitors'  filename];
        fullpath = [filepath '\' num2str(ii) '\monitors'  filename];
        fullpath = [filepath '\monitors'  filename];
        m_dot(ii) = ReadExportFile(fullpath);
    else
        filename = ['\' geo '_' fluid '-scaled_mdot=' num2str(m_dot(ii)) '_dpdl.txt'];
%         fullpath = [filepath filename];
        fullpath = [filepath '\' num2str(ii) '\solution\monitors'  filename];
        dpdl(ii) = ReadExportFile(fullpath);
    end
     
end




%% read .cas file

% filepath = 'C:\Users\alexabus\OneDrive_NTNU\SWwd\FLUENT\AdWell\Turbulence\remote\solved';
% 
% 
% for ii=1:length(dpdl)
% %     ii=4;
%     
%     filename = ['\p_' fluid '-scaled_m-dot=' num2str(m_dot(ii)) '.cas'];
%     
%     LinesToSkip = 1554;
%  
%     fid = fopen ( [filepath filename] ); % This opens a file
%     if fid ~= -1
%     for j=1:LinesToSkip 
%         fgetl ( fid ); % anytime you use fgetl you read a line and next time you use it reads next   
%     end
%     Str = fgetl ( fid );
%     fclose ( fid );
%     end
%     Key   = 'derivative';
%     Index = strfind(Str, Key);
%     dpdl(ii,1) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);    
% end

   