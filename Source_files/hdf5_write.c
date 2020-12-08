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
    // *******************code test 
    
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
    printf("malloc stuff \n"); 

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

    int info = MPI_COMM_NULL; 

    MPI_Comm_size(cartcomm, &size);

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_get(cartcomm, ndim, dims, &periods, coords); 

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
    ierr = H5Pset_fapl_mpio(plist_id, cartcomm, info);
    if(ierr<0)
    {
        printf("unsuccessful");
        exit(0); 
    }

    // create file hdf5.dat collectively 
    file_id = H5Fcreate ("hdf5.dat", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Fclose (file_id);    

    // create data space for dataset 
    // hid_t H5Screate_simple( int rank, const hsize_t * current_dims, const hsize_t * maximum_dims ); // creates a new simple dataspace and opens it for access
    H5Screate_simple(ndim, dimsf, count); // MAX dims???

    // create dataset with default properties 
    // H5Dcreate( hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t dcpl_id )
    // hid_t data_def = H5Dcreate( file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_id );   
    H5Dcreate1(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, dset_id); 
    H5Sclose(filespace); 

    // each process defines dataset in memory and writes to hyperslab in the file. 
    H5Screate_simple(ndim, count, maxdims); 
    
    //   ! Select hyperslab in the file.
    //   CALL h5dget_space_f(dset_id, filespace, ierr)
    H5Dget_space(dset_id); 
    
    //   CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, block, count, NULL); 
          
    //   ! Create property list for collective dataset write
    //   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
    plist_id = H5Pcreate(H5P_DATASET_XFER); // sets data transfer mode.

    //   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // sets data transfer mode.

    //   Write the dataset collectively. 

    // herr_t H5Dwrite( hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void * buf )
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, iodata);

    // ! Write the dataset independently. Comment out the collective
    // ! and use the following call to investigate non-collective performance.
    // !    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, arraygsize, ierr, &
    // !                     file_space_id = filespace, mem_space_id = memspace)
    
    //   Close dataspaces.
    //   CALL h5sclose_f(filespace, ierr)
    H5Sclose(filespace); 
    //   CALL h5sclose_f(memspace, ierr)
    H5Sclose(memspace); 
    
    //   Close the dataset and property list.
    H5Dclose(dset_id); 
    //   CALL h5dclose_f(dset_id, ierr)
    H5Dclose(plist_id); 
    //   CALL h5pclose_f(plist_id, ierr)
    
    //   ! Close the file.
    //   CALL h5fclose_f(file_id, ierr)
    H5Fclose(file_id);
    
    //   ! Close FORTRAN predefined datatypes.
    //   CALL h5close_f(ierr)
}