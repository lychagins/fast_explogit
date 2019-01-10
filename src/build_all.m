clear mex;

% We will place the binaries there; add this folder to the search path
addpath('../build');

% Standard exploded logit algorithm
mex('explogit.c', 'explogit_c.c', '-outdir', '../build', '-v');