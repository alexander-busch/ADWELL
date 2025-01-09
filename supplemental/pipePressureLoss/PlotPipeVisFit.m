function [  ] = PlotPipeVisFit( KPrime, nPrime, color )
%PlotPipeVisFit Summary of this function goes here
%   Detailed explanation goes here
    

NeSR = logspace(1,3);

plot(NeSR,mean(KPrime).*NeSR.^mean(nPrime),'Color',color);

clear NeSR
        
end
