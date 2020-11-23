#include <stdio.h>
#include <mpi.h>
#include "benchclock.h"
#include <stdlib.h> 

main(int argc, char **argv) 
{  
  int ierr;
  ierr = MPI_Init(&argc, &argv);
  printf("first\n");
  int firstcall = 1; // true value 
  const int numiolayer = 4; // can be constants 
  const int maxlen = 64; 
  const int numrep = 10; 
  char filedir[100]; 
  char filename[100]; 
  char iolayername[numiolayer][100]; 
  char iostring[numiolayer][100];
  int iolayer, irep; 
// Set local array size - global sizes l1, l2 and l3 are scaled
// by number of processes in each dimension
  const int n1 = 256;
  const int n2 = 256;
  const int n3 = 256;
  const int ndim = 3;
  int i1, i2, i3, j1, j2, j3, l1, l2, l3, p1, p2, p3; 
  int iodata[257][257][257]; // n1+1=257 
  int rank, size, dblesize;

//vector<int> coords(ndim); //initialise dims with 0s
  
  int dims[] = {0}; 
  int coords[ndim]; 
  const int iounit = 12; 
  const int mib = 1024*1024;
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

   // MPI_Status status;
  MPI_Comm comm = MPI_COMM_WORLD;

  int nprocs, myrank ;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
// Set 3D processor grid

  MPI_Dims_create(size, ndim, dims);

// Reverse dimensions as MPI assumes C ordering (this is not essential)?  

  p1 = dims[2];
  p2 = dims[1]; 
  p3 = dims[0];

// Compute global sizes

  l1 = p1*n1;
  l2 = p2*n2;
  l3 = p3*n3;

  MPI_Type_size(MPI_DOUBLE_PRECISION, &dblesize);  // returns number of bytes occupied by entries 
 
  mibdata = dblesize*n1*n2*n3*p1*p2*p3/mib; 

  if (rank == 0) 
  {
     printf("\n") ;
     printf("Simple Parallel IO benchmark") ;
     printf("----------------------------");
     printf("\n") ;    
     printf("Running on %c", size);
     printf("process(es)");
     printf("Process grid is ( %c", p1); printf(", %c", p2 ); printf(", %c", p3); printf(")" );
     printf("Array size is   ( %c", n1); printf(", %c", n2 ); printf(", %c", n3); printf(")" );
     printf("Global size is  ( %c", l1); printf(", %c", l2 ); printf(", %c", l3); printf(")" );
     printf("\n") ;   
     printf("Total amount of data = %c", mibdata); printf("MiB"); 
     printf("\n") ;   
     double tick = benchtick(firstcall)*1.0e6;
     printf("Clock resolution is %lf", tick); printf(", usecs");     
  }

  dims[0] = p1;
  dims[1] = p2;
  dims[2] = p3;

  MPI_Comm cartcomm; 

  MPI_Cart_create(comm, ndim, dims, periods, reorder, &cartcomm);  // new communicator to which topology information is added. 

// delete halo values and then initialise again with -1
  int i, j, k; 
  for(i = 0; i < 257; i++) {// n1 =256, n1+1 = 257
  
    for(j = 0; j < 257; j++)
    {
      for(k = 0; k < 257; k++)
      {
            iodata[i][j][k] = -1; // initialise all values with -1 
      }
    }
  }

  MPI_Cart_coords(cartcomm, rank, ndim, coords); 

  for (i3 = 0; i3 < n3; i3++)
  {
  for (i2 = 0; i2<n3; i2++)
      {
      for (i1 = 0; i1<n3; i1++) 
        {            
           j1 = coords[0]*n1 + i1;
           j2 = coords[1]*n2 + i2;
           j3 = coords[2]*n3 + i3;

           iodata[i1][i2][i3] = (j3-1)*l1*l2 + (j2-1)*l1 + j1;  
        }
      }
  }
  for (iolayer = 1; iolayer < numiolayer; iolayer++)
  { //  Skip layer if support is not compiled in
    // Expects iolayers in order: serial, MPI-IO, HDF5, NetCDF
    #ifndef WITH_SERIAL
     if (iolayer == 1) 
      { 
          continue; 
      }
    #endif
    #ifndef WITH_MPIIO
     if (iolayer == 2) 
      { 
         continue; 
      }
    #endif
    #ifndef WITH_HDF5
     if (iolayer == 3) 
      { 
          continue; 
      }
    #endif
    #ifndef WITH_NETCDF
     if (iolayer == 4) 
      { 
          continue;  
      }
    #endif
  }
  
  if(rank == 0)
  {
    printf("\n") ;
    printf(" -------"); 
    printf(iostring[iolayer-1]);  // iolayer does not start from 0
    printf(" -------");     
    printf("\n") ;
  }

  strcpy(filename, filedir); 
  strcat(filename, "/"); 
  strcat(filename, iolayername[iolayer]); 
  // string addition to get full path, trim is for trailing blankspaces

  if (rank == 0)
  {
        printf("Writing to %c", filename) ; 
        mintime = 0;
        maxiorate = 0;
        avgtime = 0;
        avgiorate = 0; 
  }
  for (irep = 1; irep<=numrep; irep++) //should it start with 0? 
  {   
    MPI_Barrier(MPI_COMM_WORLD);  
    if (rank == 0)
    {       t0 = benchtime(firstcall);   
    }
    
    switch(iolayer)
    {
      case 1:
       serialwrite(*filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined? 
      break;
      case 2:
       mpiiowrite(*filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined?
      break;
//       case 3:
//       call hdf5write(filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined?
//       break;
//       case 4:
//      call netcdfwrite(filename, iodata, n1, n2, n3, cartcomm); // function calls need to be defined? 
//      break;
      default:
      printf("Illegal value of iolayer = %d", iolayer);     
    }

    MPI_Barrier(MPI_COMM_WORLD); 
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
      printf("time = %d", time); printf(", rate = %d", iorate); printf(" MiB/s"); 
    }
  } 
  
  if (rank == 0) 
  {
    printf("mintime = %d", mintime); printf(", maxrate = %d", maxiorate); printf(" MiB/s"); 
    printf("avgtime = %d", avgtime); printf(", avgrate = %d", avgiorate); printf(" MiB/s"); 
    printf("Deleting: %c", filename);
    printf("\n"); 
  }
  
  if (rank == 0)
  {
    printf("\n"); 
    printf("--------");
    printf("Finished");
    printf("--------");    
    printf("\n");  
  }
ierr = MPI_Finalize();
return 0;

}

