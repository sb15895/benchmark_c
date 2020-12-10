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

void hdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm)
{   // *******************code test 
    
    // hid_t file_id;
    // herr_t status;
    // file_id = H5Fcreate ("file.h5", H5F_ACC_TRUNC,
    // H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Fclose (file_id);

    // *******************code test end 
    int ndim = 3; 
    // char* filename; 
    //  integer(hsize_t), dimension(ndim) :: dimsf  ! dataset dimensions.  
    char dsetname[8]; 
    strcpy(dsetname, "IntArray"); 
    // integer(hid_t) :: file_id       ! file identifier  
    int i, ierr, rank, size, periods, initialized; 
    int arraygsize[NDIM], arraystart[NDIM]; 
    int ncid, varid, oldmode, x_dimid, y_dimid, z_dimid, dimids[NDIM]; 

    // Dynamically allocate memory using malloc() 
    int* dims;
    int* coords; 
    dims = (int*)malloc(ndim * sizeof(int)); 
    coords = (int*)malloc(ndim * sizeof(int)); 

    hid_t file_id; // file identifier
    hid_t dset_id; // dataset identifier 
    hid_t filespace; // dataspace identifier in file 
    hid_t memspace; // dataspace identifier in memory
    hid_t plist_id; // property list identifier 

    hsize_t count[NDIM]; 
    hssize_t offset[NDIM]; 
    hsize_t dimsf[NDIM]; 
    hsize_t maxdims[10]; // randomly selected maximum dimension number. Not sure. 
    hsize_t block[NDIM]; // block array determines size of element blcok selected from dataspace. 
    herr_t status;

    // int info = MPI_COMM_NULL;
 
    MPI_Comm_size(cartcomm, &size);
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, ndim, dims, &periods, coords); 
    // MPI_Info info; 

    int arraysize[] = {n1+2, n2+2, n3+2}; 

    // subtract halo from array subsize 

    int arraysubsize[] = {n1, n2, n3};

    for (i = 0; i < ndim ; i++) // vector multiplication of arraysize+dims, arrayssubsize+coords
    {   
        arraygsize[i] = arraysubsize[i] * dims[i]; // arraygsize(:) = arraysubsize(:) * dims(:)
        arraystart[i] = arraysubsize[i] * coords[i];    
    }

    // copy arraygsize into dimsf 
    memcpy(dimsf, arraygsize, sizeof(int)*ndim);

    // initialise count and offset arrays 
    count[0] = n1; 
    count[1] = n2; 
    count[2] = n3; 
    
    for (i = 0; i < ndim ; i++) // defines offset used in the HDF5 file
    {   
        offset[i] = coords[i] * count[i]; // arraygsize(:) = arraysubsize(:) * dims(:)
    }

    // initialise hdf5
    
    // setup file access property list with parallel IO access
    plist_id = H5Pcreate(H5P_FILE_ACCESS); 
    H5Pset_fapl_mpio(plist_id, cartcomm, MPI_INFO_NULL);
    printf("h5pset fapl mpio \n"); 

    // create file hdf5.dat collectively
    // plist_id is fapl_id, access property list id 
    // file_id is file identifier 
    file_id = H5Fcreate ("hdf5.dat", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose (plist_id);    
    printf("file id \n"); 

    // create data space for dataset 
    filespace = H5Screate_simple(ndim, dimsf, NULL); 
    printf("H5Screate_simple \n"); 

    // dset_id = H5Pcreate(H5P_FILE_ACCESS); 
    // printf("Dataset ID \n"); 

    // create dataset with default properties 
    dset_id = H5Dcreate(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
    printf("H5Dcreate dset_id \n"); 
    H5Sclose(filespace); 

    // each process defines dataset in memory and writes to hyperslab in the file. 
    memspace = H5Screate_simple(ndim, count, NULL); // memspace is dataspace identifier
    printf ("dataspace simple \n");

    //   Select hyperslab in the file.
    filespace = H5Dget_space(dset_id); 
    printf ("dget space \n");
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL); 
    printf("Hyperslab \n"); 
          
    //  Create property list for collective dataset write
    plist_id = H5Pcreate(H5P_DATASET_XFER); // sets data transfer mode.
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // sets data transfer mode.
    printf("H5Pset dxpl mpi \n"); 

    //   Write the dataset collectively. 
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, iodata);
    printf("H5D write \n"); 
    // status = H5Dclose (dset_id);    
    // status = H5Pclose (plist_id);
    // status = H5Sclose (filespace);
    // status = H5Sclose (memspace);
    // printf("File close \n"); 
    // status = H5Fclose (file_id);
}