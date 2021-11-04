// calculate the numerical integration of \int_0^1 (4.0 / (1 + x^2))
// midpoint rule
// should approximately equals \pi
// block distribution

#include<iostream>
#include<omp.h>

#define NTHREADS 6

static long nx = 1000000000;

int main()
{
    int real_threads;
    double pi = 0.0;
    double start_time, end_time, run_time;
    double delta_x;
    double sum[NTHREADS] = {0.0};

    delta_x = 1.0 / (double) nx;

    omp_set_num_threads(NTHREADS);

    start_time = omp_get_wtime();
    #pragma omp parallel
    {
        double x;
        int id = omp_get_thread_num();
        int number_threads = omp_get_num_threads();
        int block = nx / number_threads;
        int start = id * block;
        int end = (id + 1) * block;

        if (id == number_threads - 1)
        {
            end = nx;
            real_threads = number_threads;
        }
        
        for (int i = start; i < end; i ++)
        {
            x = (i + 0.5) * delta_x;
            sum[id] += (4.0 / (1.0 + x * x));
        }
        sum[id] = delta_x * sum[id];
    } // end of parallel region
    for (int i = 0; i < NTHREADS; i++)
    {
        pi += sum[i];
    }
    end_time = omp_get_wtime();
    
    std::cout << "Pi = " << pi << std::endl;
    std::cout << "use " << real_threads << " threads, " << "take " << nx <<"steps, " 
                << end_time - start_time << "secs." << std::endl;
    
    return 0;
}


/*
use 1 threads, take 1000000000steps, 6.034secs.
use 2 threads, take 1000000000steps, 8.923secs.
use 3 threads, take 1000000000steps, 6.48secs.
use 6 threads, take 1000000000steps, 4.434secs.
use 12 threads, take 1000000000steps, 3.37secs.
*/