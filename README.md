# fast_explogit
Fast mixed exploded logit with agent-specific choice sets.

## Building the binaries ##

### Linux ###
To compile the code and save the binaries in the `build` folder, run
```
make
```
By default, the build uses Matlab's own BLAS and doesn't rely on OpenMP. To enable OpenMP support, run
```
make CFLAGS_EXT=-fopenmp LDFLAGS_EXT=-fopenmp
```
To use an external BLAS library, run
```
make BLAS_NAME=blas_name BLAS_PATH=blas_path
```
For example, assuming that OpenBLAS is installed as `~/OpenBLAS/lib/libopenblas.so`, type 
```
make BLAS_NAME=openblas BLAS_PATH=~/OpenBLAS/lib CFLAGS_EXT=-fopenmp LDFLAGS_EXT=-fopenmp
```
to enable both OpenMP and OpenBLAS.

Run 
```
make test
```
to perform optional tests. To clean up, run `make clean`.


### Windows ###
Run `src/build_all.m` from the Matlab command prompt. Run optional tests in the `tests/` folder. The build relies on Matlab's own BLAS; OpenMP support is not available.

## Usage ##

[See pdf documentation](docs/formulas.pdf)

### Example ###

Suppose that agent A chooses 2 most preferred items from a choice set {1, 3, 4, 6}. Agent B chooses 1 most preferred item from {2, 5, 6}. A's first and second best choices are 3 and 6 respectively, while B's first best is 2.

A's and B's choice sets and stated preferences determine the ordering of columns in X. Columns 1–4 represent the choice set of A, in the reverse preference order. That is, covariates of the most preferred item, 3, are placed to column 4; the second-best choice, 6, goes to column 3. Columns 1 and 2 correspond to choices 1 and 4 (their relative order is unimportant). The block in columns 5–7 correspond to choices of agent B: item 2 is in column 7, items 5 and 6 are in columns 5–6.