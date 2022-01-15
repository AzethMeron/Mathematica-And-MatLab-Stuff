close all;
clear all;
clc;

% I think there's some error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of car
l = 2;
d = 1.5;
c = 1;

% initial conditions
f0 = [ 1 2 deg2rad(45) deg2rad(15) deg2rad(0) ]; % [ x y fi theta ]

% Controls - task A,B
u = { 
    @(t) 1
    @(t) 0 
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
pause_time = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
tspan = [ tmin tmax ];
[t, sol] = ode45(@(t,f) car(t, f, l, d, u), tspan, f0, options );

% Interpolation
it = [ tmin : ts : tmax ]';
isol = interp1(t, sol, it);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Car body
Nx = isol(:,1);
Ny = isol(:,2);
Rx = Nx + d*cos(isol(:,3));
Ry = Ny + d*sin(isol(:,3));
Fx = Rx + l*cos(isol(:,4));
Fy = Ry + l*sin(isol(:,4));
% Rear wheel
rear_Px = Nx + c*sin(isol(:,3));
rear_Py = Ny - c*cos(isol(:,3));
rear_Dx = Nx - c*sin(isol(:,3));
rear_Dy = Ny + c*cos(isol(:,3));
% Middle wheel
mid_Px = Rx + c*sin(isol(:,4));
mid_Py = Ry - c*cos(isol(:,4));
mid_Dx = Rx - c*sin(isol(:,4));
mid_Dy = Ry + c*cos(isol(:,4));
% Front wheel
front_Px = Fx + c*sin(isol(:,5) + isol(:,4));
front_Py = Fy - c*cos(isol(:,5) + isol(:,4));
front_Dx = Fx - c*sin(isol(:,5) + isol(:,4));
front_Dy = Fy + c*cos(isol(:,5) + isol(:,4));

Animate( [0, 20, 0, 20],@(i,t) "Racer","x","y",0.01,it,{
    @(i,t) plot([rear_Px(i); rear_Dx(i)],[rear_Py(i); rear_Dy(i)],'-b'); 
    @(i,t) plot([mid_Px(i); mid_Dx(i)],[mid_Py(i); mid_Dy(i)],'-b'); 
    @(i,t) plot([front_Px(i); front_Dx(i)],[front_Py(i); front_Dy(i)],'-b'); 
    @(i,t) plot([Nx(i); Rx(i)],[Ny(i); Ry(i)],'-b'); 
    @(i,t) plot([Rx(i); Fx(i)],[Ry(i); Fy(i)],'-b'); 
})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = Generator(t, t0, controls, index)
    i = 1 + mod(floor(t/t0), length(controls));
    u = controls{i}(index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [df] = car(t, f, l, d, u)
    % f = [ x y alpha beta gamma ]
    df = zeros(5,1);
    df(1) = cos(f(5)) * cos(f(3)) * u{1}(t);
    df(2) = cos(f(5)) * sin(f(3)) * u{1}(t);
    df(3) = sin(f(5))*u{1}(t) / l;
    df(4) = cos(5) * sin(f(3)-f(4))*u{1}(t) / d;
    df(5) = u{2}(t);
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

