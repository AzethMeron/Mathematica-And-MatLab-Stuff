close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 2

[sol, t, delta] = solv(1,1);
delta
figure(); hold on; grid on; plot(t,sol(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 3

[sol, t] = task3(0.1);
figure(); hold on; grid on; plot(t,sol(:,1));
[sol, t] = task3(0);
figure(); hold on; grid on; plot(t,sol(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 4

[sol, t] = task4([0 1 0]);
figure(); hold on; grid on; plot(t,sol(:,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 5

[sol, t] = task5(3, 1);
figure(); hold on; grid on; plot(t,sol(:,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = odefcn2(t,y,p,q) % TASK 1
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = -p*y(2) - q*y(1);
end

function [sol, t, delta] = solv(p,q) % TASK 1
    delta = p*p - 4*q;
    tspan = [0 10];
    y0 = [1 0];
    [t, sol] = ode45(@(t,y) odefcn2(t,y,p,q), tspan, y0);
end

function dx = odefcn3(t,x,k) % TASK 2
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -x(1) - k*x(2);
end

function [sol, t] = task3(k) % TASK 2
    tspan = [0 20];
    x0 = [1 0];
    [t, sol] = ode45(@(t,x) odefcn3(t,x,k), tspan, x0);
end

function dx = odefcn4(t,x)
    dx = zeros(3,1);
    dx(1) = -4*x(1) + 2*x(2) +5*x(3);
    dx(2) = 6*x(1) - x(2) - 6*x(3);
    dx(3) = -8*x(1) + 3*x(2) + 9*x(3);
end

function [sol, t] = task4(x0)
    tspan = [0 10];
    [t, sol] = ode45(@(t,x) odefcn4(t,x), tspan, x0);
end

function xd = odefcn5(t,x,k,m)
    xd = zeros(2,1);
    xd(1) = x(2);
    xd(2) = (-k*x(1))/m;
end

function [sol, t] = task5(k, m)
    tspan = [0,10];
    x0 = [0 1];
    [t, sol] = ode45(@(t,x) odefcn5(t,x,k,m), tspan, x0);
end