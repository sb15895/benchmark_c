#include <mpi.h>
#include <string.h>
double benchtick(int a);
double benchtime(int a);
void mpiiowrite(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm); 
void serialwrite(double* iodata); 
void hdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm )
