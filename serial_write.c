#include <omp.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

//************************************SERIAL WRITE******************************************** 

void serialwrite(long int iodata[][257][257])
{
        int iounit = 10;  
        int i, j, k; 
        FILE* in_file; 
        chdir("benchio_files");   
        printf("change directory \n");   
        in_file = fopen("serial.dat", "w"); // read only  

        if (! in_file ) // equivalent to saying if ( in_file == NULL ) 
            {  
                printf("oops, file can't be read\n"); 
                exit(-1); 
            } 
        for (i = 0; i<257; i++)
            {   
                for (j = 0; j<257; j++)
                {
                    for (k =0; k<257; k++)
                    {
                        fprintf(in_file, "%d \n", iodata[i][j][k]);
                    }
                }
            } 
        printf("Data written in iodata is %d \n",iodata[200][200][200]); 
        printf("File written \n");           

        fclose(in_file); 

}
