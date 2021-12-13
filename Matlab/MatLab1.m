close all;
clear all;

% Jakub Grzana 241530

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 1
A = [1 4 ; -2 0];
B = [0 -6; -8 2];

3*A - 0.5*B;
A';
A*B;
B*A;
A*A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 2
M = [ 5 -1; 7 8 ];
det(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 3
a = [ 0:2:8;
    1:2:9;
    2:2:10;
    3:2:11;
    4:2:12];

b = [ 0:1:4;
    1:1:5;
    2:1:6;
    3:1:7;
    4:1:8;
    ];
b.*b;

c = diag(ones(1,5));
C = num2cell(c,2);
C{2};

d = [ 1 2 3 4 ];
normalised_d = d / norm(d);

e = flip(a,2);

f = a;
for i = 1:5
    f([i],[2 4]) = f([i],[4 2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 4
null = zeros(5);
unit = ones(5);
diagonal = diag([1:1:5]);
diagonal2 = diag([1:1:5].*[1:1:5]);

m = diagonal2;
m([2],[2]) = 7;

m = diagonal2;
m([1 2],:) = m([2 1],:);

m = diagonal2;
m(:, [3 5]) = m(:, [5 3]);

m = unit;
m([2], :) = 0;

m = unit;
m([1:2],[1:2]) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 5
[plus minus] = fun(10,12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TASK 6
a = input('Pwease give me fiwst vawue uwu :3\n$ ');
b = input('Pwease give me second vawue uwu :3\n$ ');
s = a+b;
fprintf("Sum of those two iws: %d\n", s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function definitions in a script must appear at the end of the file.
function [first, second] = fun(a,b)
    first = a+b;
    second = a-b
end
