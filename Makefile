RED			='\033[1;31m'
NC			='\033[0m' # No Color

CC			= gcc

CCFLAGS		= -Wall -O3 -fopenmp -v \
			-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free \
			-I$(HOME)/include -L$(HOME)/lib
			
CLDFLAGS 	= 	$(HOME)/lib/libamdlibm.a \
			$(HOME)/lib/libopenblas.a \
			-Wl,-rpath=$(HOME)/lib \
			-ltcmalloc
			
PCCFLAGS	= $(CCFLAGS) \
			-DWITHGPERFTOOLS -g

PLDFLAGS	= $(CLDFLAGS) \
			-lprofiler

# optim: explogit.c wrapper_c.c
	# $(CC) $(CCFLAGS) explogit.c wrapper_c.c	-fprofile-generate -o wrapper_c $(CLDFLAGS)
	# ./wrapper_c
	# $(CC) $(CCFLAGS) explogit.c wrapper_c.c	-fprofile-use -o wrapper_c $(CLDFLAGS)
	
# profile:
	# $(CC) $(PCCFLAGS) -c explogit.c
	# $(CC) $(PCCFLAGS) -c wrapper_c.c
	# $(CC) $(PCCFLAGS) explogit.o wrapper_c.o -o wrapper_c $(PLDFLAGS)
	# ./wrapper_c
	# pprof --text --lines ./wrapper_c explogit.profile
	
# mex: build/lcexplogit_c.o src/lcexplogit.c
	# cd src && \
	# mex CFLAGS='$$CFLAGS -fopenmp -I$(HOME)/include -std=c99' \
		# CLIBS='-L$(HOME)/lib -lgomp -lopenblas $$CLIBS' -v \
		# COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' \
		# -outdir ../build \
		# lcexplogit.c ../build/lcexplogit_c.o
		
mex:
	cd src/; \
		matlab -r "build_all; exit" -nodesktop -nosplash

test:
	cd tests/; \
		matlab -r "test_explogit; exit" -nodesktop -nojvm -nosplash

build/lcexplogit_c.o: src/lcexplogit_c.c
	$(CC) $(CCFLAGS) -c src/lcexplogit_c.c -o build/lcexplogit_c.o

clean:
	rm -f 	build/* \
		*.profile \
		*.heap \
		*.gcda
