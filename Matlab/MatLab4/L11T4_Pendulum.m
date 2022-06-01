close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of pendulum
l = [ 1 1 ];
m = [ 1 1 ];
g = 9.81;

% initial conditions
theta0 = [0 0 deg2rad(400) 0]; % [ theta1 theta2 dtheta1 dtheta2 ]
theta0_ = [0 0 deg2rad(400.1) 0]; % initial conditions for 2nd pendulum

% parameters of "simulation"
tmin = 0;
ts = 0.1;
tmax = 1000;
options = odeset('RelTol', 1e-8);
pause_time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
tspan = [ tmin tmax ];
[t, sol] = ode45(@(t,theta) double_pendulum(t, theta, l, m, g), tspan, theta0, options);
[t2, sol2] = ode45(@(t,theta) double_pendulum(t, theta, l, m, g), tspan, theta0_, options);

% Interpolation
it = [ tmin : ts : tmax ]';
isol = interp1(t, sol, it);
isol2 = interp1(t2, sol2, it);

% Coords of balls
X1 = zeros(length(it),1);
Y1 = zeros(length(it),1);
X2 = X1 + l(1)*sin(isol(:,1));
Y2 = Y1 - l(1)*cos(isol(:,1));
X3 = X2 + l(2)*sin(isol(:,2));
Y3 = Y2 - l(2)*cos(isol(:,2));

X1_ = zeros(length(it),1);
Y1_ = zeros(length(it),1);
X2_ = X1_ + l(1)*sin(isol2(:,1));
Y2_ = Y1_ - l(1)*cos(isol2(:,1));
X3_ = X2_ + l(2)*sin(isol2(:,2));
Y3_ = Y2_ - l(2)*cos(isol2(:,2));

% Animation
figure()
axes = [-3 3 -3 3];
for i = 1:length(it)
   hold off;
   plot(0,0);
   hold on;
   line = [ X1(i) Y1(i); X2(i) Y2(i) ];
   plot(line(:,1),line(:,2),'-b'); 
   line = [ X2(i) Y2(i); X3(i) Y3(i) ];
   plot(line(:,1),line(:,2),'-b'); 
   
   axis(axes);
   plot(X1(i),Y1(i),'og',X2(i),Y2(i),'og',X3(i),Y3(i),'og');
   
   % 2nd pendulum
   line = [ X1_(i) Y1_(i); X2_(i) Y2_(i) ];
   plot(line(:,1),line(:,2),'-r'); 
   line = [ X2_(i) Y2_(i); X3_(i) Y3_(i) ];
   plot(line(:,1),line(:,2),'-r'); 
   plot(X1_(i),Y1_(i),'oy',X2_(i),Y2_(i),'oy',X3_(i),Y3_(i),'oy');
   
   title(gca, sprintf('Pendulum Problem step:%d',i))
   ylabel("Y")
   xlabel("X")
   
   drawnow;
   pause(pause_time)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dtheta] = double_pendulum(t, theta, l, m, g)
    % m = [ m1, m2 ];
    % l = [ l1, l2 ];
    % g = 9.81
    % theta = [ theta1, theta2, theta1', theta2' ];
    M = [(m(1)+m(2))*l(1), m(2)*l(2)*cos(theta(1)-theta(2)); m(2)*l(1)*cos(theta(1)-theta(2)), m(2)*l(2)];
    A = [ 0, m(2)*l(2)*sin(theta(1)-theta(2)); -m(2)*l(1)*sin(theta(1)-theta(2)), 0 ];
    B = [ (m(1)+m(2))*g*sin(theta(1)); m(2)*g*sin(theta(2)) ];
    dtheta = zeros(4,1);
    dtheta(1:2) = theta(3:4);
    dtheta(3:4) = M^(-1) * (A*(theta(3:4).^2) + B) * (-1);
end