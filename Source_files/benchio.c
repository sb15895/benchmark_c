#include <stdio.h>
#include <mpi.h>
#include <stdlib.h> 
#include "bench_headerfiles.h"
#include <string.h> 
#include <hdf5.h>
#include <H5FDmpio.h> 
#include <hdf5_hl.h> 
#define N 256 
#define NDIM 3 
#define NUMIOLAYER 5

main(int argc, char **argv) 
{    
  double l1, l2, l3, p1, p2, p3, t0, t1, time, iorate, mibdata, mintime, maxiorate, avgtime, avgiorate;  
  int i1, i2, i3, j1, j2, j3, ierr, rank, dblesize, iolayer, irep, size, N1, N2, N3; 
  N1 = 256; // size of dimensions  
  N2 = 256; 
  N3 = 256;
  int firstcall = 1; // true value 
  int maxlen = 64; 
  int numrep = 1; //should be 10 ? 
  char filedir[100]; 
  char filename[100]; 
  char iolayername[NUMIOLAYER][100]; 
  char iostring[NUMIOLAYER][100];
  double* iodata;
  iodata = (double*)malloc((N+1)*(N+1)*(N+1)*sizeof(double));

  if (iodata == NULL) 
  {  
    printf("Could not allocate memory for population.\n");
    exit(0);
  }

  // Set local array size - global sizes l1, l2 and l3 are scaled
  // by number of processes in each dimension
  
  int dims[NDIM] = {0,0,0}; 
  int coords[NDIM]; 
  int iounit = 12; 
  int mib = 1024*1024;
  int reorder = 0; // false flag
  int periods[3] = {0,0,0}; //all marked false
 
  strcpy(iostring[0], "Serial");
  strcpy(iostring[1], "MPI-IO");
  strcpy(iostring[2], "HDF5");
  strcpy(iostring[3], "PHDF5");

  strcpy(iolayername[0], "serial.dat");
  strcpy(iolayername[1], "mpiio.dat");
  strcpy(iolayername[2], "hdf5.dat");
  strcpy(iolayername[3], "phdf5.dat");

  // MPI initialisation, MPI Comm, rank and process number
  ierr = MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Set 3D processor grid
  MPI_Dims_create(size, NDIM, dims);

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
     printf("Total amount of data = %f MiB \n", mibdata); 
     double tick = benchtick(firstcall)*1.0e6;
     printf("Clock resolution is %e, u secs \n", tick);  
  }

  dims[0] = p1;
  dims[1] = p2;
  dims[2] = p3;

// new communicator to which topology information is added. 
  MPI_Comm cartcomm; 
  MPI_Cart_create(comm, 3, dims, periods, reorder, &cartcomm);  
  MPI_Cart_coords(cartcomm, rank, NDIM, coords); 

// delete halo values and then initialise again with -1
  int i, j, k; 
  for(i = 0; i <=N1; i++) 
  {// n1 =256, n1+1 = 257
    for(j = 0; j <= N2; j++)
    {
      for(k = 0; k <= N3; k++)
      {
            *(iodata + i*(int)(N2*N3) + j*(int)N3 + k) = -1; // initialise all values with -1 
      }
    }
  }  

// iodata values filled in

  for(i = 0; i <= N1; i++) {
  
    for(j = 0; j <= N2; j++)
    {
      for(k = 0; k <= N3; k++)
      {
           j1 = coords[0]*N1 + i;
           j2 = coords[1]*N2 + j;
           j3 = coords[2]*N3 + k;

           *(iodata + i*(int)(N2*N3) + j*(int)N3 + k) = (j3-1)*l1*l2 + (j2-1)*l1 + j1;  
      }
    }
  }  
  i = 100; j = 100; k = 100; 
  // printf("iodata[100][100][100] after data entered %f \n", *(iodata + (int)(i*65536) + (int)(j*256) + k)); 

  // Check for directory "benchio_files"
  
  int check; 
  char* dirname = "benchmark_c/benchio_files"; 
  check = mkdir(dirname,0777); 
  if (!check) 
        printf("Directory created\n"); 
  else { 
        printf("Directory already present \n"); 
    } 
 
  // Loop to call serial write, mpiwrite, hdf5 write and adios write. 

  for (iolayer = 4; iolayer < 5; iolayer++)
  { 

    // HAVE TO COMPLETE PROCESSOR PRE DEFINITION CHECKS. 

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

    strcpy(filename, iolayername[iolayer]); 

    // reset all time parameters 
    if (rank == 0)
    {
      printf("\nWriting to benchio_files/%s \n", filename) ; 
      mintime = 0;
      maxiorate = 0;
      avgtime = 0;
      avgiorate = 0; 
    }
    for (irep = 0; irep <numrep; irep++)
    {
      // MPI_Barrier(cartcomm);  
      if (rank == 0)
      { 
              t0 = benchtime(firstcall);   
      }

      switch(iolayer)
      {
        case 0:
        if(rank == 0)
        {
          serialwrite(iodata);
          printf("Serial write completed\n\n"); 
        } 
        break;
        case 1:
          // printf("cartcomm for MPI %d \n", cartcomm); 
          mpiiowrite(iodata, N1, N2, N3, cartcomm); 
          printf("MPI write completed\n\n"); 
        break;
        case 2:
        if (rank == 0)
        {
          hdf5write(iodata, N1, N2, N3, cartcomm); 
          printf("Serial HDF5 write completed \n\n");   
        }
        case 3: 
          phdf5write(iodata, N1, N2, N3, cartcomm);
          printf("Parallel HDF5 write completed \n\n");       
        break;

        case 4:
          adioswrite(iodata, N1, N2, N3, cartcomm); 
          printf("ADIOS write completed \n\n"); 
        break; 

        default:
        printf("Illegal value of iolayer = %d \n", iolayer);     
      }
      
      // MPI_Barrier(cartcomm); 
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
        printf("Time = %e, Rate = %e MiB/s \n", time, iorate);  
        
      }
    }
  }
  
  if (rank == 0) 
  {
    printf("mintime = %e, maxrate = %e MiB/s \n", mintime, maxiorate); 
    printf("avgtime = %e, avgiorate = %e MiB/s \n", avgtime, avgiorate); 
    printf("Deleting: %s \n", filename);
    int del = remove(filename);
    
    if (!del)
      printf("The file is deleted successfully \n");
    else
      printf("The file is not deleted \n");

    printf("\n --------Finished-------- \n"); 
  }

  ierr = MPI_Finalize();
  free(iodata); 
  return 0;

}

