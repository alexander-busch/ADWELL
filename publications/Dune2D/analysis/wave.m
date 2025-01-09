
clear all;
close all;
clc;


Length = 0.4;
t_max = 10;

deltat = 0.01;
deltax = 0.01;


FREQ=0.1;
AMPL_INF=0.001;
k = 30 % Wavenumber
cp = 3 % Sigmoid centerpoint
g = 2 % Sigmoind growth

figure;

subplot(1,2,1)
x=linspace(0,Length,1/deltax);
t=linspace(0,t_max,1/deltat);
[x,time] = meshgrid(x,t);
y = AMPL_INF./(1+exp(-g.*(time-cp))) .*sin(k.*x-2.*pi.*FREQ.*time); % Compute displacement
surface(x,time,y)
view(-70,40)

subplot(1,2,2)
x = linspace(0,Length,100);
y = zeros(1,100);
time=0;
s = plot(x,y);
axis([0 Length -AMPL_INF AMPL_INF])
for i=1:t_max/deltat
    time = i*deltat % Current time
    y = AMPL_INF/(1+exp(-g*(time-cp))) * sin(k*x-2*pi*FREQ*time); % Compute displacement
    s.YData = y; % Replace plot y values
    pause(deltat) % pause to control animation speed
end;

