#include <mpi.h>
#include <string.h>
#define N 256
double benchtick(int a);
double benchtime(int a);
void mpiiowrite(long int iodata[][N+1][N+1], int n1, int n2, int n3, MPI_Comm cartcomm); 
//void serialwrite(std::string filename, int iodata[][4][4], int n1, int n2, int n3, int cartcomm); 
void serialwrite(long iodata[][N+1][N+1]); 

// #endif  // end BENCHCLOCK_HPP