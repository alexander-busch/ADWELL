
%% Clean-up
data(1,:)=[]; % Delete first row in data which were headings in Excel
condition = data(:,3) == 0; % Get rows where PV = 0
data(condition,:)=[]; % Remove respective rows
data=data(~any(isnan(data),2),:);  % Remove respective rows

%% Create unique IDs
ID = data(:,end); % Create vector containing all wellbore IDs
[Order,index] = unique(ID); % Create matrix containing unique IDs & original index
UniqueID = ID(sort(index)); % Create vector containing unique IDs with original order

%% Plot
fig = figure();
hold on;
Shearrate = logspace(-2,3);
l = length(UniqueID);
% Loop through all unique Wellbore IDs
for i = 1:l
    Wellbore = find(ID == UniqueID(i)); % Get positions of current wellbore in data matrix 
    a = Wellbore(1); % Establish lower boundary
    b = Wellbore(end); % Establish upper boundary
    A = data(a:b,:); % Generate reduced data matrix for current wellbore
%     condition = A(:,3) == 0; % Get rows where PV = 0
%     A(condition,:)=[]; % Remove respective rows
    Depth = A(:,1);
    YP = A(:,4); % Get YP
    YP = YP(:,ones(1,length(Shearrate))); % Create Yield stress matrix
    Shearstress = YP + A(:,3) ./1000 * Shearrate; % Compute ShearStress-ShearRate = f(Depth)


   
    % Plot single ShearStress-ShearRate-functions
    D = ones(length(Shearrate),1);
    for j = 1:length(Depth)
%         plot3(Shearrate,Depth(j)*D,Shearstress(j,:)./Shearrate);
        plot(Shearrate,Shearstress(j,:)./Shearrate);
    end

% surf(Shearrate,Depth,Shearstress);       
end

%    YP = data(:,4); % Get YP
%     YP = YP(:,ones(1,length(Shearrate)));
% plot(Shearrate, (YP + data(:,3)./1000 * Shearrate)./Shearrate);


set(gca,...
    'xlim', [0 1000],...
    'ylim', [1e-2 1e1],...
    'XScale','log',...
    'YScale','log',...
    'box','on',...
    'FontSize',24);

view(3) % az = –37.5, el = 30
view(2) % az = 0, el = 90.