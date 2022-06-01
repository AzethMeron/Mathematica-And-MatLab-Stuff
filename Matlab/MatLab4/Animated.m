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
tmax = 100;
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

Animate([-3 3 -3 3], @(i,t) sprintf("Sample number %d",i), "t", "y(t)", 0, it, { 
    @(i,t) plot([X1(i); X2(i)],[Y1(i); Y2(i)],'-b'); 
    @(i,t) plot([X2(i); X3(i)],[Y2(i); Y3(i)],'-b'); 
    @(i,t) plot(X1(i),Y1(i),'og',X2(i),Y2(i),'og',X3(i),Y3(i),'og');
    
    @(i,t) plot([X1_(i); X2_(i)],[Y1_(i); Y2_(i)],'-r'); 
    @(i,t) plot([X2_(i); X3_(i)],[Y2_(i); Y3_(i)],'-r'); 
    @(i,t) plot(X1_(i),Y1_(i),'og',X2_(i),Y2_(i),'og',X3_(i),Y3_(i),'og');  
})

function Animate(axes, title_, xlabel_, ylabel_, pausetime, samples, objs)
    % axes = [ x_min x_max y_min y_max ]
    % title = @(i, t) where i is sample number, and t is value of sample
    % xlabel, ylabel self explanatory
    % pausetime = number, interval between frames
    % samples = vector, your "x" (or "t")
    % objs = cell array with anonymous functions @(i,t) plot(...);
    % example: @(i,t) plot([X1(i); X2(i)],[Y1(i); Y2(i)],'-b'); 
    
    figure();
    for i = 1:length(samples)
        hold off;
        x_mid = (axes(2)+axes(1))/2;
        y_mid = (axes(4)+axes(3))/2;
        plot(x_mid, y_mid); % cleaning prev frame
        
        hold on;
        axis(axes);
        title(title_(i,samples(i)));
        xlabel(xlabel_);
        ylabel(ylabel_);
        
        for j = 1:length(objs)
            func = objs{j};
            func(i, samples(i));
        end
        
        drawnow;
        pause(pausetime); 
    end
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