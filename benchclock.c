// benchclock.cpp containing benchtick and benchtime functions
#include "benchclock.h"
#include <time.h>
#include <stdio.h>      
//using namespace std::chrono;

//chrono::high_resolution_clock Clock;

double benchtick(int firstcall)
{   double rate = 7473190834065087488; // from benchlock.f90 kind = selected int kind (18)
    double count = 7473190834065087488;
    double benchtick2; 
    double ticktime = 0.0;
    if(firstcall != 0) //if firstcall is non 0 then: 
    {
        firstcall = 0; //false
        //  call system_clock(count, rate)
           
        // ticktime = 1.0d0/dble(rate) ! what is 1.0d0?  
        ticktime = 1/rate; 
        
    }
    benchtick2 = ticktime;
    return(benchtick2); 
}


double benchtime(int firstcall)
{   double count = 7473190834065087488;
    clock_t start_t, end_t, total_t;
    start_t = clock();
    end_t = clock();
    total_t = end_t-start_t;
    double ticktime = 0.0;
    printf("Elapsed time is %ld", total_t); 
    printf("\n"); 
    if (firstcall) // if firstcall is non 0, that means its been called for the first time? 
    {
        double dummy = benchtick(firstcall); //whats the purpose?
    }


    double benchtime2 = count*ticktime; 
    return (benchtime2); 
}