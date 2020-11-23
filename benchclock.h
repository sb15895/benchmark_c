#include <mpi.h>
double benchtick(int a);
double benchtime(int a);
void mpiiowrite(char filename, int iodata[][257][257], int n1, int n2, int n3, MPI_Comm cartcomm); 
//void serialwrite(std::string filename, int iodata[][4][4], int n1, int n2, int n3, int cartcomm); 
void serialwrite(char filename, int iodata[][257][257], int n1, int n2, int n3, MPI_Comm cartcomm); 

// #endif  // end BENCHCLOCK_HPP