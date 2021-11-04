// calculate the numerical integration of \int_0^1 (4.0 / (1 + x^2))
// midpoint rule
// should approximately equals \pi
// cyclic distribution
// using critical function to avoid cache false sharing

#include<iostream>
#include<omp.h>

#define NTHREADS 18

static long nx = 2000000000;

int main()
{
    int real_threads;
    double pi = 0.0;
    double start, end, run_time;
    double delta_x;
    double sum1, sum2;
    double sum = {0.0};

    delta_x = 1.0 / (double) nx;

    omp_set_num_threads(NTHREADS);

    start = omp_get_wtime();
    #pragma omp parallel
    {
        double x;
        int id = omp_get_thread_num();
        int number_threads = omp_get_num_threads();
        double sum_par = 0.0;

        for (int i = id; i < nx; i += number_threads)
        {
            x = (i + 0.5) * delta_x;
            sum_par += (4.0 / (1.0 + x * x));
        }
        sum_par = delta_x * sum_par;

        #pragma omp critical
        {
            pi += sum_par;
        }

        if(id == 0) real_threads = number_threads;
    } // end of parallel region

    end = omp_get_wtime();
    
    std::cout << "Pi = " << pi << std::endl;
    std::cout << "use " << real_threads << " threads, " << "take " << nx <<"steps, " 
                << end - start << "secs." << std::endl;
    
    return 0;
}


/*
use 1 threads, take 1000000000steps, 6.031secs.
use 2 threads, take 1000000000steps, 3.04secs.
use 3 threads, take 1000000000steps, 2.034secs.
use 6 threads, take 1000000000steps, 1.047secs.
use 12 threads, take 1000000000steps, 0.651secs.
use 24 threads, take 1000000000steps, 0.596secs.
*/