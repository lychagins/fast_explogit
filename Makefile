RED			='\033[1;31m'
NC			='\033[0m' # No Color

CC			= gcc

CCFLAGS		= -Wall -O3 -fopenmp \
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


explogit: explogit.o wrapper_c.o
	$(CC) $(CCFLAGS) explogit.o wrapper_c.o -o wrapper_c $(CLDFLAGS)

optim: explogit.c wrapper_c.c
	$(CC) $(CCFLAGS) explogit.c wrapper_c.c	-fprofile-generate -o wrapper_c $(CLDFLAGS)
	./wrapper_c
	$(CC) $(CCFLAGS) explogit.c wrapper_c.c	-fprofile-use -o wrapper_c $(CLDFLAGS)
	
profile:
	$(CC) $(PCCFLAGS) -c explogit.c
	$(CC) $(PCCFLAGS) -c wrapper_c.c
	$(CC) $(PCCFLAGS) explogit.o wrapper_c.o -o wrapper_c $(PLDFLAGS)
	./wrapper_c
	pprof --text --lines ./wrapper_c explogit.profile
	
mex: explogit.c explogit_mex.c
	mex CFLAGS='$$CFLAGS -fopenmp -I$(HOME)/include -std=c99' \
		CLIBS='-L$(HOME)/lib -Wl,-rpath=$(HOME)/lib -lgomp -lamdlibm -lopenblas -ltcmalloc $$CLIBS' -v \
		COPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' \
		explogit_mex.c explogit.c
	echo -e ${RED}Before running Matlab, set LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6${NC}
	
wrapper_c.o: wrapper_c.c
	$(CC) $(CCFLAGS) -c wrapper_c.c

explogit.o: explogit.c
	$(CC) $(CCFLAGS) -c explogit.c

clean:
	rm -f 	explogit.o \
		wrapper_c.o \
		wrapper_c \
		explogit_mex.mexa64 \
		explogit.mexa64 \
		*.profile \
		*.heap \
		*.gcda
