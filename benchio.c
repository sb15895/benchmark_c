#include <stdio.h>
#include <mpi.h>
#include <stdlib.h> 
#include "benchclock.h"
#include <string.h> 
#define N 256 
#define NDIM 3 

long int iodata[N+1][N+1][N+1];

main(int argc, char **argv) 
{    
  int ierr; 
  double l1, l2, l3, p1, p2, p3, N1, N2, N3; 
  int i1, i2, i3, j1, j2, j3; 
  N1 = 256;
  N2 = 256;
  N3 = 256;
  printf("first\n");
  int firstcall = 1; // true value 
  int numiolayer = 4; // can be constants 
  int maxlen = 64; 
  int numrep = 10; 
  char filedir[100]; 
  char filename[100]; 
  char iolayername[numiolayer][100]; 
  char iostring[numiolayer][100];
  int iolayer, irep; 

// Set local array size - global sizes l1, l2 and l3 are scaled
// by number of processes in each dimension
  int rank, size, dblesize;

//vector<int> coords(ndim); //initialise dims with 0s
  
  int dims[NDIM] = {0,0,0}; 
  int coords[NDIM]; 
  int iounit = 12; 
  int mib = 1024*1024;
  int reorder = 0; // false flag
  int periods[3] = {0,0,0}; //all marked false
  double t0, t1, time, iorate, mibdata; 
  double mintime, maxiorate, avgtime, avgiorate; 

  strcpy(iostring[0], "Serial");
  strcpy(iostring[1], "MPI-IO");
  strcpy(iostring[2], "HDF5");
  strcpy(iostring[3], "NetCDF");

  strcpy(iolayername[0], "serial.dat");
  strcpy(iolayername[1], "mpiio.dat");
  strcpy(iolayername[2], "hdf5.dat");
  strcpy(iolayername[3], "netcdf.dat");

  strcpy(filedir, "benchio_files");

  ierr = MPI_Init(&argc, &argv);

   // MPI_Status status;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// Set 3D processor grid
  MPI_Dims_create(size, NDIM, dims);
  
  printf("hello from rank %d and process %d \n", rank, size); 
// Reverse dimensions as MPI assumes C ordering (this is not essential)?  

  p1 = dims[2];
  p2 = dims[1]; 
  p3 = dims[0];

// Compute global sizes

  l1 = p1*N1;
  l2 = p2*N2;
  l3 = p3*N3;

  MPI_Type_size(MPI_DOUBLE_PRECISION, &dblesize);  // returns number of bytes occupied by entries 
 
  mibdata = dblesize*N1*N2*N3*p1*p2*p3/mib;  

  if (rank == 0) 
  {
     printf("\n") ;
     printf("Simple Parallel IO benchmark") ;
     printf("---------------------------- \n");
     printf("Running on %d process(es) \n", size);
     printf("Process grid is ( %d, %d, %d ) \n", p1, p2, p3); 
     printf("Array size is ( %d, %d, %d ) \n", N1, N2, N3); 
     printf("Global size is ( %d, %d, %d ) \n", l1, l2, l3); 
     printf("Total amount of data = %d MiB \n", mibdata); 
     double tick = benchtick(firstcall)*1.0e6;
     printf("Clock resolution is %e, u secs \n", tick);  
  }

  dims[0] = p1;
  dims[1] = p2;
  dims[2] = p3;

  MPI_Comm cartcomm; 
  printf("cartcomm \n");

  MPI_Cart_create(comm, 3, dims, periods, reorder, &cartcomm);  // new communicator to which topology information is added. 

  printf("cart create \n"); 
// delete halo values and then initialise again with -1
  int i, j, k; 
  for(i = 0; i <= N; i++) {// n1 =256, n1+1 = 257
  
    for(j = 0; j <= N; j++)
    {
      for(k = 0; k <= N; k++)
      {
            iodata[i][j][k] = -1; // initialise all values with -1 
      }
    }
  }
  printf("iodata initialised \n"); 

  MPI_Cart_coords(cartcomm, rank, NDIM, coords); 

  printf("cart coords  \n"); 

  for (i3 = 0; i3 <= N3; i3++)
  {
  for (i2 = 0; i2 <= N2; i2++)
      {
      for (i1 = 0; i1 <= N1; i1++) 
        {            
           j1 = coords[0]*N1 + i1;
           j2 = coords[1]*N2 + i2;
           j3 = coords[2]*N3 + i3;

           iodata[i1][i2][i3] = (j3-1)*l1*l2 + (j2-1)*l1 + j1;  
        }
      }
  }

  printf("iodata filling \n"); 

  // for (iolayer = 0; iolayer < numiolayer; iolayer++)
  // { //  Skip layer if support is not compiled in
  //   // Expects iolayers in order: serial, MPI-IO, HDF5, NetCDF
  //   #ifndef WITH_SERIAL
  //    if (iolayer == 1) 
  //     { 
  //         continue; 
  //     }
  //   #endif
  //   #ifndef WITH_MPIIO
  //    if (iolayer == 2) 
  //     { 
  //        continue; 
  //     }
  //   #endif
  //   #ifndef WITH_HDF5
  //    if (iolayer == 3) 
  //     { 
  //         continue; 
  //     }
  //   #endif
  //   #ifndef WITH_NETCDF
  //    if (iolayer == 4) 
  //     { 
  //         continue;  
  //     }
  //   #endif
  // }
  
  // if(rank == 0)
  // {
  //   printf("\n") ;
  //   printf(" -------"); 
  //   printf(iostring[iolayer-1]);  // iolayer does not start from 0
  //   printf(" -------");     
  //   printf("\n") ;
  // }

  iolayer = 0; 
  strcpy(filename, filedir); 
  strcat(filename, "/"); 
  strcat(filename, iolayername[iolayer]); 
  // string addition to get full path, trim is for trailing blankspaces

  printf("filename \n"); 

  if (rank == 0)
  {
        printf("Writing to %s \n", filename) ; 
        mintime = 0;
        maxiorate = 0;
        avgtime = 0;
        avgiorate = 0; 
  }
  // for (irep = 1; irep<=numrep; irep++) //should it start with 0? 
  // for (irep = 1; irep<=1; irep++) //should it start with 0? 
  // {   
  //   MPI_Barrier(MPI_COMM_WORLD);  
  //   if (rank == 0)
  //   {       t0 = benchtime(firstcall);   
  //   }
    
    // if (rank == 0) 
    // {
    serialwrite(*filename, iodata);
    // } 
    printf("Serial write completed\n");
//     switch(iolayer)
//     {
//       case 1:
//        serialwrite(*filename, iodata, N1, N2, N3, cartcomm); // function calls need to be defined? 
//       break;
//       case 2:
//        mpiiowrite(*filename, iodata, N1, N2, N3, cartcomm); // function calls need to be defined?
//       break;
// //       case 3:
// //       call hdf5write(filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined?
// //       break;
// //       case 4:
// //      call netcdfwrite(filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined? 
// //      break;
//       default:
//       printf("Illegal value of iolayer = %d", iolayer);     
//     }

    // MPI_Barrier(MPI_COMM_WORLD); 
    if(rank == 0)
    {
    double t1 = benchtime(firstcall); 
    time = t1-t0; 
    iorate = mibdata/time;
    avgtime = avgtime + time/numrep;
    avgiorate = avgiorate + iorate/numrep;

      if (maxiorate < iorate)
      {
        maxiorate = iorate;
        mintime = time;
      }
      printf("time = %e, rate = %e MiB/s \n", time, iorate);  
    }
  // } 
  
  if (rank == 0) 
  {
    printf("mintime = %e, maxrate = %e MiB/s \n", mintime, maxiorate); 
    printf("avgtime = %e, avgiorate = %e MiB/s \n", avgtime, avgiorate); 
    printf("Deleting: %s", filename);
    printf("\n"); 
    printf("--------");
    printf("Finished");
    printf("--------");    
    printf("\n");  
  }

ierr = MPI_Finalize();
return 0;

}

