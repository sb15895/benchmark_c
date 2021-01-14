#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <adios2_c.h>
#define ENGINE "BPFile"

void adioswrite(double* iodata, int N1, int N2, int N3, MPI_Comm cartcomm)
{   int rank, size; 
    MPI_Comm_rank(cartcomm, &rank);
    MPI_Comm_size(cartcomm, &size);
    
    // Init the ADIOS subsystem
    adios2_error adiosErr;
    
    // adios2_adios* adios = adios2_init(MPI_COMM_WORLD, adios2_debug_mode_on);
    adios2_adios* adios = adios2_init(MPI_COMM_WORLD);
    adios2_io* ioWrite = adios2_declare_io(adios, "A");
    adiosErr = adios2_set_engine(ioWrite, ENGINE);
    const char greeting[] = "Hello World from ADIOS2";
    
    // writer(adios, greeting);
    adios2_io* io = adios2_declare_io(adios, "hello-world-writer");
    adios2_variable* var_greeting = adios2_define_variable(io, "Greeting", adios2_type_string, 0, NULL, NULL, NULL, adios2_constant_dims_true);
    adios2_engine* engine = adios2_open(io, "hello-world-c.bp", adios2_mode_write);
    adios2_put(engine, var_greeting, greeting, adios2_mode_deferred);
    adios2_close(engine);
    adios2_finalize(adios);
}
