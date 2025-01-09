function [ fig ] = CreateFigure( fig_title, label_x, label_y, axis, format)
%	CreateFigure Create and format figure 
%   Detailed explanation goes here

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
                            
fig = figure('Name',fig_title,'color','w','Units','centimeters','Position',format);

hold on;
grid on;
box on;

% title(fig_title);
xlabel(label_x); ylabel(label_y); 

set(gca,'XScale',axis{1},'YScale',axis{2});

end

