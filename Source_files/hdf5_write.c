#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <hdf5.h>
#define N 256 
#define NDIM 3 

void hdf5write(double* iodata, int n1, int n2, int n3, MPI_Comm cartcomm)
{
    // code test 
    hid_t file_id;
    herr_t status;
    file_id = H5Fcreate ("file.h5", H5F_ACC_TRUNC,
    H5P_DEFAULT, H5P_DEFAULT);
    status = H5Fclose (file_id);

    // int ndim = 3; 
    // char* filename; 
    // //  integer(hsize_t), dimension(ndim) :: dimsf  ! dataset dimensions.  
    // char dsetname[8]; 
    // strcpy(dsetname, "IntArray"); 
    // // integer(hid_t) :: file_id       ! file identifier  
    // int i, ierr, rank, size, periods; 
    // int dims[ndim], coords[ndim], arraygsize[ndim], arraystart[ndim]; 
    // int ncid, varid, oldmode, x_dimid, y_dimid, z_dimid, dimids[ndim]; 

    // hid_t file_id; 
    // hid_t dset_id;
    // hid_t filespace; 
    // hid_t memspace; 
    // hid_t plist_id; 
    // hsize_t count[ndim]; 
    // hssize_t offset[ndim]; 
    // herr_t status;

    // MPI_Comm_size(cartcomm, &size);
    // MPI_Comm_rank(cartcomm, &rank);

    // MPI_Cart_get(cartcomm, ndim, dims, &periods, coords); 

    // int arraysize[] = {n1+2, n2+2, n3+2}; 

    // // subtract halo from array subsize 

    // int arraysubsize[] = {n1, n2, n3};

    // for (i = 0; i < ndim ; i++) // vector multiplication of arraysize+dims, arrayssubsize+coords
    // {   
    //     arraygsize[i] = arraysubsize[i] * dims[i]; // arraygsize(:) = arraysubsize(:) * dims(:)
    //     arraystart[i] = arraysubsize[i] * coords[i];    
    //     // printf("rank %d i %d arraygsize %d arraystart %d \n", myrank, i, arraygsize[i], arraystart[i]);
    // }

    // int dimsf[ndim]; 

    // // copy arraygsize into dimsf 
    // memcpy(dimsf, arraygsize, sizeof(int)*ndim);

    // // initialise count and offset arrays 
    // count[0] = n1; 
    // count[1] = n2; 
    // count[2] = n3; 
    
    // for (i = 0; i < ndim ; i++) // defines offset used in the HDF5 file
    // {   
    //     offset[i] = coords[i] * count[i]; // arraygsize(:) = arraysubsize(:) * dims(:)
    //     // printf("rank %d i %d arraygsize %d arraystart %d \n", myrank, i, arraygsize[i], arraystart[i]);
    // }

    // // initialise hdf5
    
    // // setup file access property list with parallel IO access?
    // // CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    // // CALL h5pset_fapl_mpio_f(plist_id, cartcomm, info, ierr)

    // // create file hdf5.dat collectively 
    // chdir("benchio_files"); 
    // file_id = H5Fcreate ("hdf5.dat", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Fclose (file_id);    

    // // create data space for dataset 


}