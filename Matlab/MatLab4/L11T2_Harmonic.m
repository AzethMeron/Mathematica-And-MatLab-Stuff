close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of oscillator
m = 1;
dzeta = [0.5 1 1.5 2];
k = 0.5;

% initial conditions
x0 = [0 1];

% parameters of "simulation"
tmin = 0;
ts = 0.01;
tmax = 30;
options = odeset('RelTol', 1e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% You can just calculate and plot in the same loop.
% I'm dividing it for clarity

results = {}; % vessek for results
it = [ tmin : ts : tmax ]';

for i = 1:length(dzeta)
   % solving
    tspan = [ tmin tmax ];
    c = dzeta(i) * 2 * sqrt(m*k);
    [t, sol] = ode45(@(t,x) func(t,x, k, m, c), tspan, x0, options);

    % interpolating
    isol = interp1(t, sol, it);
    
    % keeping result
    results{i} = isol(:,1) ;
end

figure();
hold on; grid on;
xlabel('t');
ylabel('x(t)');
for i = 1:length(dzeta)
    plot(it, results{i}, 'DisplayName', sprintf('Dzeta = %.2f',dzeta(i)));
end
lgd = legend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dt] = func(t, x, k, m, c)
    % x = [ x(t), x'(t) ]
    dt = zeros(2,1);
    dt(1) = x(2);
    dt(2) = (-c*x(2) - k*x(1)) / m;
end