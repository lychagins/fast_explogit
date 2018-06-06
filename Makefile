CC	= gcc

# CCFLAGS	= -Wall -O3 -fopenmp \
	  # -L$(HOME)/lib -I$(HOME)/include \
	  # -Wl,-rpath=$(HOME)/lib
CCFLAGS	= -Wall -O3 -fopenmp \
	-DWITHGPERFTOOLS -g -L$(HOME)/lib -I$(HOME)/include \
	-Wl,-rpath=$(HOME)/lib \
	-Wl,--no-as-needed -lprofiler -ltcmalloc -Wl,--as-needed 

explogit: explogit.o wrapper_c.o
	$(CC) $(CCFLAGS) explogit.o wrapper_c.o -o wrapper_c -lm

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
