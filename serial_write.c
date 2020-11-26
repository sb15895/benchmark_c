#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

//************************************SERIAL WRITE******************************************** 

void serialwrite(double* iodata)
{
    int iounit = 10;  
    int i, j, k; 
    int N1, N2, N3; 
    N1 = 256; 
    N2 = 256; 
    N3 = 256;
    FILE* in_file; 
    chdir("benchio_files");   
    printf("Directory changed \n");   
    in_file = fopen("serial.dat", "w"); // read only  
    printf("File opened \n");
    if (! in_file ) // equivalent to saying if ( in_file == NULL ) 
    {  
        printf("oops, file can't be read\n"); 
        exit(-1); 
    } 
    i = 100; j = 100; k = 100; 
    printf("Value of iodata[100][100][100] after passed to function %f \n", *(iodata + (int)(i*65536) + (int)(j*256) + k)); 


    for (i = 0; i < N1; i++) 
    {   
        for (j = 0; j < N2; j++)
        {
            for (k = 0; k < N3; k++)
            {
                // fprintf(in_file, "%d \n", *(iodata + i*(int)(N2*N3) + j*(int)(N3) + k));
                fprintf(in_file, "%f ", *(iodata + (int)(i*65536) + (int)(j*256) + k));
            }
        }
    } 
    printf("File written \n");           
    fclose(in_file); 

}
