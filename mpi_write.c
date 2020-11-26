#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include<string.h>
#define N 256 
#define NDIM 3 

// void mpiio()
// {
// #ifdef WITH_MPIIO 
// use mpi  
// #endif 
// }
                       
//************************************MPI WRITE******************************************** 

void mpiiowrite(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm)
{   
    // printf("function called \n"); 
    int ndim = 3;
    int i, j, k, initialized, finalized, arraygsize[ndim], arraystart[ndim], ierr, size, comm, rank, periods, argc, nprocs, myrank; 
    MPI_File fh; 
    MPI_Status status;
    int info = 0;      
    int order = 0; 

    //means to initialise different dimesions with 0s till n1+1, etc
    for (i = 0; i<n1+1; i++)
    {  
        for (j = 0; j<n2+1; j++)
        {
            for(k = 0; k<n3+1; k++)
            { 
                *(iodata + (int)(i*65536) + (int)(j*256) + k) = 0; 
            }
        }
    }

    // printf("iodata filled \n"); 

    // Dynamically allocate memory using malloc() 
    int* dims;
    int* coords; 
    dims = (int*)malloc(ndim * sizeof(int)); 
    coords = (int*)malloc(ndim * sizeof(int)); 

    // printf("malloc stuff \n"); 
  
    // Check if the memory has been successfully 
    // allocated by malloc or not 
    if (dims == NULL || coords == NULL) { 
        printf("Memory not allocated.\n"); 
        exit(0); 
    } 

    // MPI should be initialised only once. 

    MPI_Initialized(&initialized);
    if (!initialized)
    {
    MPI_Init(NULL, NULL);
    printf("MPI initialised \n");
    }

    MPI_Comm_size(cartcomm, &nprocs);
    MPI_Comm_rank(cartcomm, &myrank);
    // printf("MPI rank %d MPI processor %d\n", myrank, nprocs);
    
    MPI_Cart_get(cartcomm, ndim, dims, &periods, coords); 
  
    int arraysize[] = {n1+2, n2+2, n3+2};  

    // Halo subtract from array subsize
    int arraysubsize[] = {n1, n2, n3}; 

    for (i = 0; i < ndim ; i++) // vector multiplication of arraysize+dims, arrayssubsize+coords
    {   
        arraygsize[i] = arraysubsize[i] * dims[i]; 
        arraystart[i] = arraysubsize[i] * coords[i];    
        printf("rank %d i %d arraygsize %d arraystart %d \n", myrank, i, arraygsize[i], arraystart[i]);
    }

    // MPI_Type_create_subarray creates an MPI dataype to extract data and write it to file. 
    MPI_Datatype mpi_subarray, filetype; 
    //define filetype for process: portion of global array that process owns.
    MPI_Type_create_subarray(ndim, arraygsize, arraysubsize, arraystart, MPI_ORDER_C, MPI_DOUBLE, &filetype); // order?     
    MPI_Type_commit(&filetype);
    // printf("first subarray \n");     
    
    //define subarray for process: portion of local array thats to be written. 
    for (i = 0; i < ndim; i++)
    {
        arraystart[i] = 1; 
    }
    MPI_Type_create_subarray(ndim, arraysize, arraysubsize, arraystart, MPI_ORDER_C, MPI_DOUBLE, &mpi_subarray); // order?     
    MPI_Type_commit(&mpi_subarray); 
    // printf("second subarray \n");     

    //write files to directory benchio_files
    chdir("benchio_files"); 
    // printf("Directory changed \n"); 
    
    // open file "mpiio.dat" using write only mode and create if the file doesnt exist.
    MPI_File_open(cartcomm, "mpiio.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh); 
    printf("MPI file open \n"); 

    // set view using filetype data type and write file
    MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    printf("MPI set view \n"); 
    //remove halo data by passing MPI subarray type
    MPI_File_write_all(fh, iodata, 1, mpi_subarray, &status);
    printf("MPI file write all \n"); 

    MPI_File_close(&fh);  // close file 
    printf("MPI close file \n"); 
    
    // MPI_Finalized(&finalized);
    // if (!finalized)
    MPI_Finalize();

    free(dims);
    free(coords); 
}
                           

