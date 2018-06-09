CC			= gcc

CCFLAGS		= -Wall -O3 -fopenmp \
			-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free \
			-I$(HOME)/include -L$(HOME)/lib
			
LDFLAGS 	= 	$(HOME)/lib/libamdlibm.a \
			$(HOME)/lib/libopenblas.a \
			-Wl,-rpath=$(HOME)/lib \
			-ltcmalloc
			
PCCFLAGS	= $(CCFLAGS) \
			-DWITHGPERFTOOLS -g

PLDFLAGS	= $(LDFLAGS) \
			-lprofiler

explogit: explogit.o wrapper_c.o
	$(CC) $(CCFLAGS) explogit.o wrapper_c.o	-o wrapper_c $(LDFLAGS)

profile:
	$(CC) $(PCCFLAGS) -c explogit.c
	$(CC) $(PCCFLAGS) -c wrapper_c.c
	$(CC) $(PCCFLAGS) explogit.o wrapper_c.o -o wrapper_c $(PLDFLAGS)
	./wrapper_c
	pprof --text --lines ./wrapper_c explogit.profile
	
mex: explogit.c explogit_mex.c
	mex -v CFLAGS="$(CFLAGS) -fopenmp" -lgomp explogit_mex.c explogit.c

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
		*.heap
