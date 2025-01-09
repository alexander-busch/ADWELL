function [c_D] = cD_Morrison2016( Re_p )
%	cD_Morrison2013 Drag law of Morrison(2013)
%   Detailed explanation goes here

% Drag coefficient
c_D = (24./Re_p)+(2.6.*(Re_p./5))./(1+(Re_p./5).^1.52)+(0.411.*(Re_p./263000).^(-7.94))./(1+(Re_p./263000).^(-8))+0.25.*(Re_p./1e6)./(1+Re_p./1e6);

end % of function