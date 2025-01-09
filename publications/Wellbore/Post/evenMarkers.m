function [xmarkers,ymarkers,h] = evenMarkers(x,y,NumMarkers)
  xmarkers = logspace(log10(x(1)),log10(x(end)),NumMarkers);
  ymarkers = interp1(x,y,xmarkers);
  h = plot(xmarkers,ymarkers,'ko'); 
end