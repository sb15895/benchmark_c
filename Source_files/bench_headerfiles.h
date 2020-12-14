#include <mpi.h>
#include <string.h>
double benchtick(int a);
double benchtime();
void serialwrite(double* iodata); 
void mpiiowrite(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm); 
void hdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm); 
void phdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm); 