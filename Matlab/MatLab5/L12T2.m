close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of car
l = 2;
d = 1;

% initial conditions
f0 = [ 1 2 deg2rad(45) deg2rad(0) ]; % [ x y fi theta ]


% Controls - task A,B
u = { 
    @(t) 4 
    @(t) 0.1 
};


%{
% Controls - task C
t0 = 0.2;
controls = {
  [1;0]
  [0;1]
  [1;0]
  [0;-1]
};
u = { 
    @(t) Generator(t, t0, controls, 1)
    @(t) Generator(t, t0, controls, 2)
};
%}

%{
% Controls - task D
u = { 
    @(t) 0.5*sin(2*t + 5) 
    @(t) 0.7*sin(1*t + 15) 
};
%}

% parameters of "simulation"
tmin = 0;
ts = 0.1;
tmax = 100;
options = odeset('RelTol', 1e-5);
pause_time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
tspan = [ tmin tmax ];
[t, sol] = ode45(@(t,f) car(t, f, l, u), tspan, f0, options );

% Interpolation
it = [ tmin : ts : tmax ]';
isol = interp1(t, sol, it);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Car body
rear_X = isol(:,1);
rear_Y = isol(:,2);
front_X = rear_X + l*cos(isol(:,3));
front_Y = rear_Y + l*sin(isol(:,3));
% Rear wheel
rear_wheel_LX = rear_X - d*sin(isol(:,3));
rear_wheel_LY = rear_Y + d*cos(isol(:,3));
rear_wheel_PX = rear_X + d*sin(isol(:,3));
rear_wheel_PY = rear_Y - d*cos(isol(:,3));
% Front wheel
front_wheel_LX = front_X - d*sin(isol(:,3) + isol(:,4));
front_wheel_LY = front_Y + d*cos(isol(:,3) + isol(:,4));
front_wheel_PX = front_X + d*sin(isol(:,3) + isol(:,4));
front_wheel_PY = front_Y - d*cos(isol(:,3) + isol(:,4));

Animate( [-10, 20, -10, 20], @(i,t) "Racer", "x", "y", 0.01, it, {
    @(i,t) plot([rear_wheel_LX(i); rear_wheel_PX(i)],[rear_wheel_LY(i); rear_wheel_PY(i)],'-b'); 
    @(i,t) plot([rear_X(i); front_X(i)],[rear_Y(i); front_Y(i)],'-b'); 
    @(i,t) plot([front_wheel_LX(i); front_wheel_PX(i)],[front_wheel_LY(i); front_wheel_PY(i)],'-b'); 
})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = Generator(t, t0, controls, index)
    i = 1 + mod(floor(abs(t)/t0), length(controls));
    u = controls{i}(index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [df] = car(t, f, l, u)
    % f = [ x y fi theta ]
    df = zeros(4,1);
    df(1) = cos(f(4)) * cos(f(3)) * u{1}(t);
    df(2) = cos(f(4)) * sin(f(3)) * u{1}(t);
    df(3) = sin(f(4))*u{1}(t) / l;
    df(4) = u{2}(t);
end

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

