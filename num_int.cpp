#include<iostream>
#include<omp.h>
static long nx = 1000000000;
int main()
{
    double pi = 0.0;
    double delta_x, sum, x;
    double start, end;

    omp_set_num_threads(8);

    delta_x = 1.0 / (double) nx;

    start = omp_get_wtime();
    #pragma omp parallel for default(shared) private(x)
        for (int i = 0; i < nx; i++)
        {
            x = (i + 0.5) * delta_x;
            pi+= (4.0 / (1.0 + x * x));
        }

    pi = delta_x * pi;
    end = omp_get_wtime();
    
    std::cout << "Pi = " << pi << std::endl;
    std::cout << "time = " << end - start << " secs." << std::endl;

    return 0;
}
