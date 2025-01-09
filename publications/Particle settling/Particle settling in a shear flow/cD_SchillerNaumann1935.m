function [c_D] = cD_SchillerNaumann1935( Re_p )
%cD_Rp_Renaud2004 Power law drag law of Schiller & Naumann (1935)
%   Detailed explanation goes here

% Drag coefficient
c_D = (24./Re_p).*(1+0.15.*(Re_p.^0.687));

% High Re correction
c_D(Re_p > 1000)=0.44;

end % of function