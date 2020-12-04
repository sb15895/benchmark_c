// benchclock.cpp containing benchtick and benchtime functions
#include "bench_headerfiles.h"
#include <time.h>
#include <stdio.h>     
#include <time.h>
#include <sys/time.h>

double benchtick(int firstcall)
{   double rate = 7473190834065087488; // from benchlock.f90 kind = selected int kind (18)
    double count = 7473190834065087488;
    double benchtick2; 
    double ticktime = 0.0;
    if(firstcall != 0) //if firstcall is non 0 then: 
    {
        firstcall = 0; //false
        ticktime = 1/rate;         
    }
    benchtick2 = ticktime;
    return(benchtick2); 
}


double benchtime(int firstcall)

{   
    struct timeval current_time;
    gettimeofday(&current_time, NULL);
    double timeval = current_time.tv_usec; 
    return (timeval); 
}