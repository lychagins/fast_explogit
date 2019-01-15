clear;

% Our compiled MEX binaries live here
addpath('../build');

%--------------------------------------------------------------------------
% Basic tests
%--------------------------------------------------------------------------
% Simple logit, one observation
logl_expected = 1 - log(1+exp(1));
logl_actual = lcexplogit(1, 1, [0 1], uint16(1), uint16(1));
assert(abs(logl_expected - logl_actual) < 1e-15, ...
    'Likelihood does not match U/(1-exp(U)) under simple logit.');

% Latent class exploded logit with one class has to match plain exploded
% logit.
% Simulated data, general case. 
a = cast(csvread('resources/sparse_structure.csv'), 'uint16');
nskipped = a(:,1);
nlisted = a(:,2);
a = csvread('resources/choice_x.csv');
choice_idx = cast(a(:, 1), 'uint16');
x = a(:, 2:3);

[logl, grad] = lcexplogit([1 1], 1, x', nskipped, nlisted);
assert(abs(logl + 401.4539) < 1e-4, ...
    'Loglikelihood does not match expected value on the test dataset.');
assert(all(abs(grad + [57.4651;4.5635]) < 1e-4), ...
    'Loglikelihood gradient does not match expected value on the test dataset.');

%--------------------------------------------------------------------------
% Testing multiple types case.
%--------------------------------------------------------------------------
% Three identical types. Same likelihood function as before
[logl3, grad3] = lcexplogit([2 5 1 1 1 1 1 1], 3, x', nskipped, nlisted);
assert(abs(logl3 + 401.4539) < 1e-4, ...
    'Loglikelihood does not match expected value on the test dataset.');
% Zero derivatives w.r.t. type shares
assert(all(abs(grad3(1:2)) < 1e-8), ...
    'Loglikelihood gradient does not match expected value on the test dataset.');
% Sum of derivatives w.r.t. to coefficients on X = respective coefficient
% from the plain exploded logit
assert(all(abs(sum(grad3(3:2:7)) - grad(1)) < 1e-8), ...
    'Loglikelihood gradient does not match expected value on the test dataset.');
assert(all(abs(sum(grad3(4:2:8)) - grad(2)) < 1e-8), ...
    'Loglikelihood gradient does not match expected value on the test dataset.');

disp('Latent class exploded logit: PASSED');