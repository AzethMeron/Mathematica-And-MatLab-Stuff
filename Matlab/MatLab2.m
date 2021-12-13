close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 1
p = [1 1 -6];
roots(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 2
f(5);
g(1,2);
[first second] = h(f(3), 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 3 
sine(3/2 * pi);

% TASK 4
rommel(3/2*pi);

% TASK 5
mysing(-20);

% TASK 6
mysum([1 2 3]);

% TASK 7
sumPositive([2 3 4]);

% TASK 8
[positive, negative] = sumAll([1 2 3 -5 -6 -7]);

% TASK 9
mypoly(10, [3 2 1]);

% TASK 10
mygcd(10,5);

% TASK 11
t = 0:0.1:10;
figure(); hold on; grid on; 
plot(t, arrayfun(@rommel, t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 2 functions

function output = f(x)
    output = x^3 / (x^4 + 1);
end

function output = g(x,y)
    output = sqrt(25 - x^2 - y^2);
end

function [first, second] = h(x, y)
    first = x+2*y;
    second = x*y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 3 functions

function nomine = sine(x)
    nomine = abs(sin(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 4 functions

function zhukov = rommel(x)
    smigly = sin(x);
    if smigly >= 0
        zhukov = smigly;
    else
        zhukov = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 5 functions

function sign = mysing(x)
    if x > 0
        sign = 1;
    elseif x == 0
        sign = 0;
    else
        sign = -1;
    end
    % ss = sign(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 6 functions

function ss = mysum(x)
    ss = 0;
    for i=1:length(x)
        ss = ss + x(i);
    end
    % ss = sum(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 7 functions
function xd = sumPositive(x)
    xd = 0;
    for i=1:length(x)
        if x(i) > 0
            xd = xd + x(i);
        end
    end
    % xd = sum(x(x>0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 8 functions
function [ara, uwu] = sumAll(x)
    ara = sum(x(x>0));
    uwu = sum(x(x<0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 9 functions
function val = mypoly(x, c)
    val = 0;
    for i=1:length(c)
        val = val + c(i)* x^(i-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 10 functions
function output = mygcd(m,n)
    output = gcd(m,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
