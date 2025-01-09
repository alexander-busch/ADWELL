%% Analysis of fit functions
% Rheometric data obtained at SINTEF Petroleaum AS, Bergen, 09.07.15
% (C) by Alexander Busch, NTNU, Dec. 2015

clc;
clear;


%% Fit functions
K =      0.1475;
n =      0.8942;
PL = K.*Up_SR.^(n-1);


K =      0.1475
n =      0.8941
ny_0 =   4.895e-08 
E = ny_0*ones(length(Up_SR),1) + K.*Up_SR.^(n-1);


K =     0.02592;
n =     0.5798;
ny_0 =  0.1992;
C = ny_0./(ones(length(Up_SR),1)+(K.*Up_SR).^n);


K =     0.02624;
a =     0.5811;
n =     0.4225;
ny_0 =  0.1992;
CY = ny_0./((ones(length(Up_SR),1)+(K.*Up_SR).^a)).^((1-n)/a);



figure()
hold on;

plot(Up_SR,PL,'g*')
plot(Up_SR,E,'g*')
plot(Up_SR,C,'g*')
plot(Up_SR,CY,'g*')


