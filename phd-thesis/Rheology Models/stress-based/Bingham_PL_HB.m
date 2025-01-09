clc;
clear;
close all;

%% Data

load Ezekiel2012.mat;

ShearRate = Ezekiel2012(1,:);
Fluid1 = Ezekiel2012(2,:);
Fluid2 = Ezekiel2012(3,:);
Fluid3 = Ezekiel2012(4,:);
Fluid4 = Ezekiel2012(5,:);

%% Constitutive equation/Rheology model
Bingham = 'Tau_0+K*x';
PL = 'K*x^n';
HB = 'Tau_0+K*x^n';
% Cross = 

%% Fit

% B
Results_B = zeros(4,3);

Lower = [0 0];
StartPoint = [0.913375856139019 0.63235924622541];
Upper = [100 100];

for i = 1:4
    a = coeffvalues(createFit(ShearRate,Ezekiel2012(i+1,:),Bingham, Lower, StartPoint, Upper));
    Results_B(i,1) = a(1);
    Results_B(i,2) = a(2);
end

% PL
Results_PL = zeros(4,3);

Lower = [0 0];
StartPoint = [0.913375856139019 0.0975404049994095];
Upper = [100 1];

for i = 1:4
    a = coeffvalues(createFit(ShearRate,Ezekiel2012(i+1,:),PL, Lower, StartPoint, Upper));
    Results_PL(i,1) = a(1);
    Results_PL(i,3) = a(2);
end

% HB
Results_HB = zeros(4,3);

Lower = [0 0 0];
StartPoint = [0.913375856139019 0.63235924622541 0.0975404049994095];
Upper = [100 100 1];

for i = 1:4
    a = coeffvalues(createFit(ShearRate,Ezekiel2012(i+1,:),HB, Lower, StartPoint, Upper));
    Results_HB(i,:) = a;
end


%% Plot data

figure('name','Ezekiel (2012) - Shear stress vs. Shear rate');
hold on;
for i = 1:4
    plot(Ezekiel2012(1,:), Ezekiel2012(i+1,:),'*');
end

figure('name','Ezekiel (2012) - App. viscosity vs. Shear rate');
hold on;
for i = 1:4
    plot(ShearRate, Ezekiel2012(i+1,:)./ShearRate,'*');
end
set(gca,'XScale','log','YScale','log');

%% Plot fits


for i = 1:4
B = Results_B(i,2)+Results_B(i,1).*ShearRate;
PL = Results_PL(i,1).*ShearRate.^Results_PL(i,3);
HB = Results_HB(i,2)+Results_HB(i,1).*ShearRate.^Results_HB(i,3);


figure ('name',['Fluid' num2str(i)]);

subplot(2,1,1); 
plot(Ezekiel2012(1,:), Ezekiel2012(i+1,:),'*');
plot(ShearRate, B);
plot(ShearRate, PL);
plot(ShearRate, HB);
legend ('Data', 'Bingham', 'Power-law', 'Herschel-Bulkley');

subplot(2,1,2);
hold on;
plot(Ezekiel2012(1,:), Ezekiel2012(i+1,:)./ShearRate,'*');
plot(ShearRate, B./ShearRate);
plot(ShearRate, PL./ShearRate);
plot(ShearRate, HB./ShearRate);
set(gca,'XScale','log','YScale','log');
legend ('Data', 'Bingham', 'Power-law', 'Herschel-Bulkley');

end