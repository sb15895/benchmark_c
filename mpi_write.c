#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include<string.h>
// void mpiio()
// {
// #ifdef WITH_MPIIO 
// use mpi  
// #endif 
// }
                       
//************************************MPI WRITE******************************************** 

void mpiiowrite(long int iodata[][257][257], int n1, int n2, int n3, MPI_Comm cartcomm)
{   printf("function called \n"); 
    int ndim = 3; 
    int i, j, k, initialized, finalized, arraygsize[ndim], arraystart[ndim], ierr, size, comm, mpi_subarray, rank, periods, argc; 
    MPI_File fh; 
    int info= 0; // when does it get defined? 
    //integer (kind=MPI_OFFSET_KIND) :: disp = 0
    //integer, dimension(MPI_STATUS_SIZE) :: status
    // double precision, dimension(0:n1+1,0:n2+1,0:n3+1) :: iodata // is it 0 to n1+1? 
    //means to initialise different dimesions with 0s till n1+1, etc
    for (i = 0; i<n1+1; i++)
    {   for (j = 0; j<n2+1; j++)
        {
            for(k = 0; k<n3+1; k++)
            { 
                iodata[i][j][k] = 0; 
            }
        }
    }
    printf("iodata filled \n"); 

    int* dims;
    int* coords; 

    // Dynamically allocate memory using malloc() 
    dims = (int*)malloc(ndim * sizeof(int)); 
    coords = (int*)malloc(ndim * sizeof(int)); 

    printf("malloc stuff \n"); 
  
    // Check if the memory has been successfully 
    // allocated by malloc or not 
    if (dims == NULL || coords == NULL) { 
        printf("Memory not allocated.\n"); 
        exit(0); 
    } 

    // MPI should be initialised only once. 

    MPI_Initialized(&initialized);
    if (!initialized)
    MPI_Init(NULL, NULL);

    printf("MPI initialised \n");

    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    printf("MPI rank %d \n", myrank);
    // Set 3D processor grid
    MPI_Comm_size(comm, &size);
    // dims should be 0; signals that dims create should fill in values. 
    MPI_Dims_create(size, ndim, dims);
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, ndim, dims, &periods, coords); 
    printf("cartcomm \n");

    int arraysize[] = {n1+2, n2+2, n3+2};  //n1 etc is size or values inserted? 

    // Define filetype for this process, ie what portion of the global array
    // this process owns; starting positions use C-indexing (ie counting from 0).

    int arraysubsize[] = {n1, n2, n3}; 
    // Halo subtract from array subsize

    //int arraygsize[ndim]; 
    
    for (i = 0; i < ndim ; i++) // vector multiplication of arraysize+dims, arrayssubsize+coords
    {    arraygsize[i] = arraysubsize[i] * dims[i]; 
         arraystart[i] = arraysubsize[i] * coords[i];    
    }
    int order; // not sure 

    int filetype; 
    MPI_Type_create_subarray(ndim, arraygsize, arraysubsize, arraystart, order, MPI_DOUBLE, &filetype); // order? 
    
    MPI_Type_commit(&filetype); 

    for (i =0; i<ndim;i++)
    {
        arraystart[i] = 1; 
    }

    chdir("benchio_files"); 

    MPI_Type_create_subarray(ndim, arraysize, arraysubsize, arraystart, order, MPI_DOUBLE, &mpi_subarray);
    MPI_Type_commit(&mpi_subarray);
    MPI_File_open(cartcomm, "mpiio.dat", MPI_MODE_WRONLY+MPI_MODE_CREATE, info, &fh); //not sure what mpi_mode etc etc means
    MPI_File_close(&fh);  // close file 

    MPI_Finalized(&finalized);
    if (!finalized)
    MPI_Finalize();
}
                           

