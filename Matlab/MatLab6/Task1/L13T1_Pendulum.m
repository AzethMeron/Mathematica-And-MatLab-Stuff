close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of pendulum
l = [ 1 1 ];
m = [ 1 1 ];
g = 9.81;

% initial conditions
theta0 = [ 0 0 deg2rad(400) 0]; % [ theta1 theta2 dtheta1 dtheta2 ]
theta0_ = [ 0 0 deg2rad(400.1) 0]; % [ theta1 theta2 dtheta1 dtheta2 ]

% parameters of "simulation"
tmin = 0;
ts = 0.1;
tmax = 10;
pause_time = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
sim('Pendulum');

% Interpolation
it = [ tmin : ts : tmax ]';
theta1 = interp1(tout, theta1, it);
theta2 = interp1(tout, theta2, it);
tout = it;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coords of balls
X1 = zeros(length(tout),1);
Y1 = zeros(length(tout),1);
X2 = X1 + l(1)*sin(theta1);
Y2 = Y1 - l(1)*cos(theta1);
X3 = X2 + l(2)*sin(theta2);
Y3 = Y2 - l(2)*cos(theta2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
theta0 = theta0_;
sim('Pendulum');

% Interpolation
it = [ tmin : ts : tmax ]';
theta1 = interp1(tout, theta1, it);
theta2 = interp1(tout, theta2, it);
tout = it;

% Coords of balls
X1_ = zeros(length(tout),1);
Y1_ = zeros(length(tout),1);
X2_ = X1_ + l(1)*sin(theta1);
Y2_ = Y1_ - l(1)*cos(theta1);
X3_ = X2_ + l(2)*sin(theta2);
Y3_ = Y2_ - l(2)*cos(theta2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Animation
Animate([-3 3 -3 3], @(i,t) sprintf("Double pendulum, sample %d",i), "x", "y", pause_time, tout, {
    @(i,t) plot([X1(i); X2(i)],[Y1(i); Y2(i)],'-b'); 
    @(i,t) plot([X2(i); X3(i)],[Y2(i); Y3(i)],'-b'); 
    @(i,t) plot(X1(i),Y1(i),'og', X2(i),Y2(i),'og', X3(i),Y3(i),'og');
    @(i,t) plot([X1_(i); X2_(i)],[Y1_(i); Y2_(i)],'-r'); 
    @(i,t) plot([X2_(i); X3_(i)],[Y2_(i); Y3_(i)],'-r'); 
    @(i,t) plot(X1_(i),Y1_(i),'og', X2_(i),Y2_(i),'og', X3_(i),Y3_(i),'og');
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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