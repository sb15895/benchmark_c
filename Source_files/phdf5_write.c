#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <hdf5.h>
#include <H5FDmpio.h> 
#define N 256 
#define NDIM 3 
#define FILENAME "phdf5.dat"

void phdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm)
{   char dsetname[8]; 
    strcpy(dsetname, "IntArray"); 
    int i, ierr, rank, size, periods, initialized; 
    int arraygsize[NDIM], arraystart[NDIM]; 

    // Dynamically allocate memory using malloc() 
    int* dims;
    int* coords; 
    dims = (int*)malloc(NDIM* sizeof(int)); 
    coords = (int*)malloc(NDIM* sizeof(int)); 

    hid_t file_id; // file identifier
    hid_t dset_id; // dataset identifier 
    hid_t filespace; // dataspace identifier in file 
    hid_t memspace; // dataspace identifier in memory
    hid_t plist_id; // property list identifier 

    hsize_t count[NDIM]; 
    hssize_t offset[NDIM]; 
    hsize_t dimsf[NDIM] = {n1, n2, n3};  // specifies the dimensions of dataset, dimsf[0] number of rows, dimsf[1] number of columns, dimsf[2] so on..
    hsize_t maxdims[NDIM]; // to specify maximum dimensions. 
    hsize_t block[NDIM]; // block array determines size of element blcok selected from dataspace. 
    herr_t status;

    // Obtain MPI keys 
    MPI_Comm_size(cartcomm, &size);
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, NDIM, dims, &periods, coords); 
    MPI_Info info  = MPI_INFO_NULL; 
 
    //Each process defines dataset in memory and writes it to the hyperslab in the file.  

    // initialise count and offset arrays 
    count[0] = dimsf[0]/size; // rows are divided between each process due to C memory definition. Other dimensions stay same. 
    count[1] = dimsf[1]; 
    count[2] = dimsf[2]; 

    offset[0] =rank * count[0];
    offset[1] = 0;
    offset[2] = 0; 

    // HDF5 initialisation
    // setup file access property list with parallel IO access
    plist_id = H5Pcreate(H5P_FILE_ACCESS); 
    printf("plist id create \n"); 

    H5Pset_fapl_mpio(plist_id, cartcomm, info);
    printf("fapl_mpio \n"); 

    // return file id with link to dat file "phdf5.dat"
    file_id = H5Fcreate (FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose (plist_id);    
    printf("file id \n"); 

    // return dataspace with 3 dimensions and size of each of those dimensions.
    filespace = H5Screate_simple(NDIM, dimsf, NULL); 
    printf("H5Screate_simple \n"); 

    // create dataset with default properties 
    dset_id = H5Dcreate(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
    printf("H5Dcreate dset_id \n"); 
    H5Sclose(filespace);

    // each process defines dataset in memory and writes to hyperslab in the file. 
    memspace = H5Screate_simple(NDIM, count, NULL); 
    printf ("dataspace simple \n");     
    filespace = H5Dget_space(dset_id); // makes a copy of the dataspace FOR the dataset dset_id
    printf ("dget space \n");

    //   Select hyperslab in the file.
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);      
    printf("Hyperslab \n"); 
          
    //  Create property list for collective dataset write
    plist_id = H5Pcreate(H5P_DATASET_XFER); // sets data transfer mode.
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // sets data transfer mode.
    printf("H5Pset dxpl mpi \n");   

    //   Write the dataset collectively. 
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, iodata);
    printf("H5D write \n"); 

    // Free resources

    // close dataspaces 
    status = H5Sclose (memspace);
    status = H5Sclose (filespace);

    // close dataset and property list  
    status = H5Dclose (dset_id);    
    status = H5Pclose (plist_id); 

    // close file 
    printf("File close \n"); 
    status = H5Fclose (file_id);

    // free pointers
    free(dims);
    free(coords); 
}