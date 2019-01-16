ifdef BLAS_NAME
 BLAS_LFLAGS := -L$(BLAS_PATH) -Wl,-rpath=$(BLAS_PATH) -l$(BLAS_NAME)
 BLAS_DFLAG := -DEXTERNAL_BLAS
else
 BLAS_LFLAGS := -lmwblas
endif

export BLAS_LFLAGS
export BLAS_DFLAG

ifdef LDFLAGS_EXT
export LDFLAGS_EXT
endif

ifdef CFLAGS_EXT
export CFLAGS_EXT
endif

VPATH = src:tests

source = build_all.m lcexplogit.c lcexplogit.h lcexplogit_c.c explogit.c explogit.h explogit_c.c

mex = build/explogit.mexa64 build/lcexplogit.mexa64

test = test_explogit.m test_lcexplogit.m



$(mex): $(source)
ifdef BLAS_NAME
	@echo "Building with an external BLAS"
	@echo "BLAS name: $$BLAS_NAME"
	@echo "BLAS library folder: $$BLAS_PATH"
else
	@echo "Building with the internal Matlab BLAS"
endif
	@cd src/; matlab -r "build_all; exit" -nodesktop -nosplash;

test: $(mex)
	cd tests/; matlab -r "test_explogit; test_lcexplogit; exit" -nodesktop -nojvm -nosplash

.PHONY: clean
clean:
	-rm -f $(mex)