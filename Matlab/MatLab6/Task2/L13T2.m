close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of car
l = 2;
d = 1;
c = 0.5;

% initial conditions
f0 = [ 1 2 deg2rad(45) deg2rad(0) deg2rad(0) ]; % [ x y phi1 phi0 theta ]

% Controls - task A,B
u = [2 0.1];

% parameters of "simulation"
tmin = 0;
ts = 0.05;
tmax = 30;
pause_time = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
sim('car')

% Interpolation
it = [ tmin : ts : tmax ]';
X = interp1(tout, X, it);
Y = interp1(tout, Y, it);
phi1 = interp1(tout, phi1, it);
phi0 = interp1(tout, phi0, it);
theta = interp1(tout, theta, it);
tout = it;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Car body
Rx = X;
Ry = Y;
Nx = Rx - d*cos(phi1);
Ny = Ry - d*sin(phi1);
Fx = Rx + l*cos(phi0);
Fy = Ry + l*sin(phi0);
% Rear wheel
rear_Px = Nx + c*sin(phi1);
rear_Py = Ny - c*cos(phi1);
rear_Dx = Nx - c*sin(phi1);
rear_Dy = Ny + c*cos(phi1);
% Middle wheel
mid_Px = Rx + c*sin(phi0);
mid_Py = Ry - c*cos(phi0);
mid_Dx = Rx - c*sin(phi0);
mid_Dy = Ry + c*cos(phi0);
% Front wheel
front_Px = Fx + c*sin(theta + phi0);
front_Py = Fy - c*cos(theta + phi0);
front_Dx = Fx - c*sin(theta + phi0);
front_Dy = Fy + c*cos(theta + phi0);

Animate( [0, 20, 0, 20],@(i,t) "Racer","x","y",0.01,tout,{
    @(i,t) plot([rear_Px(i); rear_Dx(i)],[rear_Py(i); rear_Dy(i)],'-b'); 
    @(i,t) plot([mid_Px(i); mid_Dx(i)],[mid_Py(i); mid_Dy(i)],'-b'); 
    @(i,t) plot([front_Px(i); front_Dx(i)],[front_Py(i); front_Dy(i)],'-b'); 
    @(i,t) plot([Nx(i); Rx(i)],[Ny(i); Ry(i)],'-b'); 
    @(i,t) plot([Rx(i); Fx(i)],[Ry(i); Fy(i)],'-b'); 
})


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

