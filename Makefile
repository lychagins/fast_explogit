CC			= gcc
CCFLAGS		= -Wall -O3 -fopenmp -I$(HOME)/include -L$(HOME)/lib
PROFFLAGS	= $(CCFLAGS) \
	-DWITHGPERFTOOLS -g  \
	-Wl,-rpath=$(HOME)/lib \
	-Wl,--no-as-needed -lprofiler -ltcmalloc -Wl,--as-needed 

explogit: explogit.o wrapper_c.o
	$(CC) $(CCFLAGS) explogit.o wrapper_c.o $(HOME)/lib/libopenblas.a -o wrapper_c -lm

profile:
	$(CC) $(PROFFLAGS) -c explogit.c
	$(CC) $(PROFFLAGS) -c wrapper_c.c
	$(CC) $(PROFFLAGS) explogit.o wrapper_c.o -o wrapper_c -lm -lopenblas
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
