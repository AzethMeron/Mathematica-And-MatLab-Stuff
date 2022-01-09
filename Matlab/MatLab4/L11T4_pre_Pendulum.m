close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of pendulum
l = 1;
m = 1;
g = 9.81;

% initial conditions
theta0 = [1 0];

% parameters of "simulation"
tmin = 0;
ts = 0.01;
tmax = 10;
options = odeset('RelTol', 1e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
tspan = [ tmin tmax ];
[t, sol] = ode45(@(t,theta) single_pendulum(t, theta, l, m, g), tspan, theta0, options);

% Interpolation
it = [ tmin : ts : tmax ]';
isol = interp1(t, sol, it);

% Coords of balls
X1 = zeros(length(it),1);
Y1 = zeros(length(it),1);
X2 = X1 + l*sin(isol(:,1));
Y2 = Y1 - l*cos(isol(:,1));

% Animation
figure()
axes = [-3 3 -3 3];
for i = 1:length(it)
   hold off;
   line = [ X1(i) Y1(i); X2(i) Y2(i) ];
   plot(line(:,1),line(:,2),'-b');
   
   hold on;
   axis(axes);
   plot(X1(i),Y1(i),'og',X2(i),Y2(i),'og');
   
   title(gca, sprintf('Pendulum Problem step:%d',i))
   ylabel("Y")
   xlabel("X")
   
   drawnow;
   pause(0.01) % pauza
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dtheta] = single_pendulum(t, theta, l, m, g)
    dtheta = zeros(2, 1);
    dtheta(1) = theta(2);
    dtheta(2) = -g/l*theta(1);
end