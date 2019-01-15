clear mex;

% We will place the binaries there; add this folder to the search path
addpath('../build');

% Standard exploded logit algorithm
mex('explogit.c', 'explogit_c.c', '-outdir', '../build', '-v');

mex COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' ...
    CFLAGS='\$CFLAGS -std=c99 \$BLAS_DFLAG' ...
    CLIBS='\$CLIBS \$BLAS_LFLAGS' ...
    -outdir ../build -v lcexplogit.c lcexplogit_c.c