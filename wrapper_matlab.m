% Read data
% beta = csvread('beta.txt');
% x = csvread('x.txt');
% nskipped = cast(csvread('nskipped.txt'), 'int32');
% nlisted = cast(csvread('nlisted.txt'), 'int32');
% save('data.mat');

clear;
load('../Data/data.mat');

tic;
[logl, grad] = explogit_mex(beta, 20, x, nskipped, nlisted);
toc;

disp(logl);