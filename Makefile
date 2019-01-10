CC			= gcc

CCFLAGS		= -Wall -O3 -v \
			-I$(HOME)/include -L$(HOME)/lib
			
CLDFLAGS 	= $(HOME)/lib/libopenblas.a \
			-Wl,-rpath=$(HOME)/lib
			
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
	
#mex: build/lcexplogit_c.o src/lcexplogit.c
	# cd src && \
	# mex CFLAGS='$$CFLAGS -fopenmp -I$(HOME)/include -std=c99' \
		# CLIBS='-L$(HOME)/lib -lgomp -lopenblas $$CLIBS' -v \
		# COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' \
		# -outdir ../build \
		# lcexplogit.c ../build/lcexplogit_c.o
		
mex:
	$(CC) $(CCFLAGS) -c src/lcexplogit_c.c -o build/lcexplogit_c.o
	mex CFLAGS='$$CFLAGS -I$(HOME)/include -std=c99' \
		CLIBS='-L$(HOME)/lib -Wl,-rpath=$(HOME)/lib -lopenblas $$CLIBS' \
		COPTIMFLAGS='-O3 -DNDEBUG' \
		LDOPTIMFLAGS='-O3' \
		-outdir build -v \
		src/lcexplogit.c build/lcexplogit_c.o
	cd src/; \
		matlab -r "build_all; exit" -nodesktop -nosplash; \
		

test:
	cd tests/; \
		matlab -r "test_explogit; test_lcexplogit; exit" -nodesktop -nojvm -nosplash

clean:
	rm -f 	build/* \
		*.profile \
		*.heap \
		*.gcda
