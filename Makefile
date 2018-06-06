CC	= gcc

CFLAGS	= -Wall \
	  -DWITHGPERFTOOLS -g -L$(HOME)/lib -I$(HOME)/include \
	  -Wl,-rpath=$(HOME)/lib \
	  -Wl,--no-as-needed -lprofiler -ltcmalloc -Wl,--as-needed 

explogit: explogit.o wrapper_c.o
	$(CC) $(CFLAGS) explogit.o wrapper_c.o -o wrapper_c -lm

mex: explogit.c explogit_mex.c
	mex explogit_mex.c explogit.c

wrapper_c.o: wrapper_c.c
	$(CC) $(CFLAGS) -c wrapper_c.c

explogit.o: explogit.c
	$(CC) $(CFLAGS) -c explogit.c

clean:
	rm -f 	explogit.o \
		wrapper_c.o \
		wrapper_c \
		explogit_mex.mexa64 \
		explogit.mexa64 \
		*.profile \
		*.heap
