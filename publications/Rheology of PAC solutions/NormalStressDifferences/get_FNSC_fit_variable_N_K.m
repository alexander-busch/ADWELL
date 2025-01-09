function [ FNSC_fit, K_coeff, n_coeff ] = get_FNSC_fit_variable_N_K( AngularFrequency,FirstNormalStressDiff,ShearRate )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




coefficients = zeros(length(AngularFrequency)-1,3);
for i=1:length(AngularFrequency)-1
    fitresult = createFit_PowerLaw(AngularFrequency(i:i+1),FirstNormalStressDiff(i:i+1));
    close gcf;
    coefficients(i,1)=(AngularFrequency(i)+AngularFrequency(i+1))/2;
    coefficients(i,2)=fitresult.K;
    coefficients(i,3)=fitresult.n;
end
clear fitresult;


% figure
% hold on
% yyaxis left
% plot(coefficients(:,1),coefficients(:,2))
% yyaxis right
% plot(coefficients(:,1),coefficients(:,3))


addpath 'E:\OneDrive_NTNU\SWwd\MATLAB\AdWell\Rheology\Models\other';
K_coeff = createfit_exp(coefficients(:,1), coefficients(:,2));
n_coeff = createfit_exp(coefficients(:,1), coefficients(:,3));
clear coefficients;
    
FNSC_fit = K_coeff.a1.*exp(-((ShearRate-K_coeff.b1)./K_coeff.c1).^2).*ShearRate.^(n_coeff.a1.*exp(-((ShearRate-n_coeff.b1)./n_coeff.c1).^2)-2);


% figure
% hold on
% yyaxis left
% plot(AngularFrequency_PAC4,K_coeff.a1.*exp(-((AngularFrequency_PAC4-K_coeff.b1)./K_coeff.c1).^2))
% yyaxis right
% plot(AngularFrequency_PAC4,n_coeff.a1.*exp(-((AngularFrequency_PAC4-n_coeff.b1)./n_coeff.c1).^2))






end

