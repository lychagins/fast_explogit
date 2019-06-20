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