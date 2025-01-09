function [ Re_cr_Bi ] = CompCrRePV( He, alliis )
%CompCrRePV Compute critical Reynolds number element-wise

Re_cr_Bi = zeros(size(He));
for ii = 1:alliis
    for jj=1:min(size(He))
        for kk=1:min(size(He))
            % Alternative: Hanks (1963) as given in Darby and Melson (1981) or
            % Chhabra and Richardson (2008)
            if He(jj,kk,ii)==0
                Re_cr_Bi(jj,kk,ii)=2100;
            elseif He(jj,kk,ii)<1700
                Re_cr_Bi(jj,kk,ii)=2100./(1+8.3e-8*log10(He(jj,kk,ii)).^13);
            elseif He(jj,kk,ii)<5e4
                Re_cr_Bi(jj,kk,ii)=80.*He(jj,kk,ii).^0.4;
            else
                Re_cr_Bi(jj,kk,ii)=25.*He(jj,kk,ii).^0.5;
            end
        end
    end
end
    
end

