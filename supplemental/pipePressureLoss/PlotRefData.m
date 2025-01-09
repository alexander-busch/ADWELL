function [  ] = PlotRefData( x, y, color )
%PlotRefData Summary of this function goes here
%   Detailed explanation goes here

plot(x, y,'--','Color',color);
plot(x, 0.9.*y,':','Color',color);
plot(x, 1.1.*y,':','Color',color);

end

