#include <mpi.h>
#include <string.h>
#define N 256
double benchtick(int a);
double benchtime(int a);
void mpiiowrite(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm); 
//void serialwrite(std::string filename, int iodata[][4][4], int n1, int n2, int n3, int cartcomm); 
void serialwrite(double* iodata); 

// #endif  // end BENCHCLOCK_HPP