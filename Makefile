# link the output files to create the executable
serial: benchio.o mpi.o benchclock.o 
	mpicc $^ -o $@

%.o: %.c
	mpicc -I. -g -o $@ -c $<

clean:
	rm -rf *.o ./serial
