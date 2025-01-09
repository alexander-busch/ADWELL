function [ fig ] = CreateFigure( fig_title, label_x, label_y, axisdefinition, format)
%	CreateFigure Create and format figure 
%   Detailed explanation goes here

if ischar(format)    

    switch format
        case 'DINA5'
            format = [1 1 21.0 14.8];
        case 'DINA4'
            format = [1 1 29.7 21.0];
        case 'DINA3'
            format = [1 1 42.0 29.7];
        case 'DINA2'
            format = [1 1 59.4 42.0];
        case 'DINA1'
            format = [1 1 84.1 59.4];
        case 'DINA0'
            format = [1 1 118.9 84.1];
        case 'PowerPoint'
            format = [1 1 26.67 14.28];
        otherwise
            
    end
    
end
                            
fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',format);

hold on;
grid on;
box on;

% title(fig_title);
xlabel(label_x); ylabel(label_y); 

if strcmp(axisdefinition,'off')
    axis off;
else
    set(gca,'XScale',axisdefinition{1},'YScale',axisdefinition{2});
end

end

