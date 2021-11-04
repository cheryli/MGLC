// calculate the numerical integration of \int_0^1 (4.0 / (1 + x^2))
// midpoint rule
// should approximately equals \pi
// cyclic distribution
// use private variable to avoid cache false sharing

#include<iostream>
#include<omp.h>

#define NTHREADS 12

static long nx = 1000000000;

int main()
{
    int real_threads;
    double pi = 0.0;
    double start, end, run_time;
    double delta_x;
    double sum[NTHREADS] = {0.0};

    delta_x = 1.0 / (double) nx;

    omp_set_num_threads(NTHREADS);

    start = omp_get_wtime();
    #pragma omp parallel
    {
        double x;
        int id = omp_get_thread_num();
        int number_threads = omp_get_num_threads();
        double sum_par = 0.0;

        // #pragma omp critical
        // std::cout << id << " " << std::hex << &sum_par << std::endl;
        /* sum_1 and sum_2  are seperated by exactly the size of 4 L1 cache: 256KB */

        for (int i = id; i < nx; i += number_threads)
        {
            x = (i + 0.5) * delta_x;
            sum_par += (4.0 / (1.0 + x * x));
        }
        sum_par = delta_x * sum_par;

        sum[id] = sum_par;

        if(id == 0) real_threads = number_threads;
    } // end of parallel region
    for (int i = 0; i < NTHREADS; i++)
    {
        pi += sum[i];
    }
    end = omp_get_wtime();
    
    std::cout << "Pi = " << pi << std::endl;
    std::cout << "use " << real_threads << " threads, " << "take " << nx <<"steps, " 
                << end - start << "secs." << std::endl;
    
    return 0;
}


/*
use 1 threads, take 1000000000steps, 6.078secs.
use 2 threads, take 1000000000steps, 3.03secs.
use 3 threads, take 1000000000steps, 2.027secs.
use 6 threads, take 1000000000steps, 1.038secs.
use 12 threads, take 1000000000steps, 0.648secs.
*/