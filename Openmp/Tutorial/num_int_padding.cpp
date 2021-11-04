// calculate the numerical integration of \int_0^1 (4.0 / (1 + x^2))
// midpoint rule
// should approximately equals \pi
// cyclic distribution
// use array padding to avoid cache false sharing

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
    double sum[NTHREADS][8] = {0.0};

    delta_x = 1.0 / (double) nx;

    omp_set_num_threads(NTHREADS);

    start = omp_get_wtime();
    #pragma omp parallel
    {
        double x;
        int id = omp_get_thread_num();
        int number_threads = omp_get_num_threads();

        for (int i = id; i < nx; i += number_threads)
        {
            x = (i + 0.5) * delta_x;
            sum[id][0] += (4.0 / (1.0 + x * x));
        }
        sum[id][0] = delta_x * sum[id][0];
        // std::cout << "test " << std::endl;
        if(id == 0) real_threads = number_threads;
    } // end of parallel region
    for (int i = 0; i < NTHREADS; i++)
    {
        pi += sum[i][0];
    }
    end = omp_get_wtime();
    
    std::cout << "Pi = " << pi << std::endl;
    std::cout << "use " << real_threads << " threads, " << "take " << nx <<"steps, " 
                << end - start << "secs." << std::endl;
    
    return 0;
}


/*
use 1 threads, take 1000000000steps, 6.038secs.
use 2 threads, take 1000000000steps, 3.048secs.
use 3 threads, take 1000000000steps, 2.062secs.
use 6 threads, take 1000000000steps, 1.084secs.
use 12 threads, take 1000000000steps, 0.685secs.
*/