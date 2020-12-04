# define c compiler 
CC = mpicc

# define compile time flags 
CFLAGS = -I/opt/ohpc/pub/libs/intel/hdf5/1.10.4/include

# define headerfiles
DEPS = bench_headerfiles.h

# define library paths in addition to /usr/lib 
LDFLAGS = -L/opt/ohpc/pub/libs/intel/hdf5/1.10.4/lib -lhdf5_hl -lhdf5 -Wl,-rpath -Wl,/opt/ohpc/pub/libs/intel/hdf5/1.10.4/lib

# define directories 
INC = -I/opt/ohpc/pub/libs/intel/hdf5/1.10.4

# define C source files 
SRCS = benchio.c benchclock.c serial_write.c mpi_write.c hdf5_write.c 

#define C object files, all .c files converted to .o files 
OBJS = $(SRCS:.c=.o)

#define executable file 
MAIN = bench

.PHONY: depend clean

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LDFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
