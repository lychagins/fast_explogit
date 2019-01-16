% Windows: Run this script to compile MEX binaries
% Linux: Use make instead of running this script directly

clear mex;

% We will place the binaries there; add this folder to the search path
addpath('../build');

% Standard exploded logit algorithm
mex('explogit.c', 'explogit_c.c', '-outdir', '../build', '-v');


if ispc
    mex COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' ...
        CFLAGS='\$CFLAGS -std=c99' ...
        -lmwblas -outdir ../build -v lcexplogit.c lcexplogit_c.c
else
    mex COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' ...
        CFLAGS='\$CFLAGS -std=c99 \$BLAS_DFLAG \$CFLAGS_EXT' ...
        CLIBS='\$CLIBS \$BLAS_LFLAGS \$LDFLAGS_EXT' ...
        -outdir ../build -v lcexplogit.c lcexplogit_c.c
end