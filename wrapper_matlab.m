% Read data
% beta = csvread('beta.txt');
% x = csvread('x.txt');
% nskipped = cast(csvread('nskipped.txt'), 'int32');
% nlisted = cast(csvread('nlisted.txt'), 'int32');
% save('data.mat');

clear;
load('../Data/data.mat');
z = x';
b = beta;
b = reshape(b(20:end), [], 20);
u = z*b;
u(1) = 213;

tic;
u = z*b;
[logl, grad] = explogit_mex(beta, 20, x, nskipped, nlisted, u);
toc;

disp(logl);