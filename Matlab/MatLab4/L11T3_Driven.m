close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of oscillator
m = 1;
c = 0.2;
k = 0.5;
omega_zero_multiplier = [ 0.5 1 1.5 ];
Fo = 1;

% initial conditions
x0 = [0 0];

% parameters of "simulation"
tmin = -50;
ts = 0.01;
tmax = 200;
options = odeset('RelTol', 1e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step response
tspan = [ tmin tmax ];
[t, sol] = ode45(@(t,x) unit_func(t,x, k, m, c), tspan, x0, options);

% Interpolating
it = [ tmin : ts : tmax ]';
isol = interp1(t, sol, it);

% Plotting
figure(); hold on; grid on;
xlabel('t');
ylabel('x(t)');
title('Damped oscillator with heaviside driving force');
plot(it, isol(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sine

results = {}; % vessek for results
it = [ tmin : ts : tmax ]';
omega_zero = sqrt(k/m);

for i = 1:length(omega_zero_multiplier)
    % solving
    tspan = [ tmin tmax ];
    omega = omega_zero_multiplier(i) * omega_zero;
    [t, sol] = ode45(@(t,x) sine_func(t,x, k, m, c, Fo, omega), tspan, x0, options);

    % interpolating
    isol = interp1(t, sol, it);
    
    % keeping result
    results{i} = isol(:,1) ;
end

figure();
hold on; grid on;
xlabel('t');
ylabel('x(t)');
for i = 1:length(omega_zero_multiplier)
    omega = omega_zero_multiplier(i) * omega_zero;
    plot(it, results{i}, 'DisplayName', sprintf('omega = %.2f * w0',omega_zero_multiplier(i)));
end
lgd = legend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dt] = unit_func(t, x, k, m, c)
    dt = zeros(2,1);
    dt(1) = x(2);
    dt(2) = (-c*x(2) - k*x(1) + heaviside(t)) / m;
end

function [dt] = sine_func(t, x, k, m, c, Fo, omega)
    dt = zeros(2,1);
    dt(1) = x(2);
    dt(2) = (-c*x(2) - k*x(1) + Fo*sin(omega*t)) / m;
end