/* This code simulates a 2-dimensional laplace equation
 * using jacobi interation
 * accelated by openacc
 * author: btli(2021)
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<omp.h>

#define index(x, y, nx) ((x) * (nx) + (y))    // a transform of 2-d array index

void initial(double *__restrict A, double *__restrict  A_new, int nx, int ny);
double jacobi(double *__restrict A, double *__restrict A_new, int nx, int ny);
void swap(double *__restrict A, double *__restrict A_new, int nx, int ny);

int main()
{
    const int nx = 4096;
    const int ny = 4096;
    const int itc_max = 1000;

    const double tolerance = 1e-6;
    double error = 1.0;

    double *restrict A = (double *)malloc(sizeof(double) * nx * ny);
    double *restrict A_new = (double *)malloc(sizeof(double) * nx * ny);

    initial(A, A_new, nx, ny);

    printf("Solve 2-D Laplace equation using jacobi relaxation\nmesh size: %d x %d\n", nx, ny);

    double start_time = omp_get_wtime();

    // jacobi interation
    int itc = 0;
    
    #pragma acc data copy(A[:nx*ny]) copyin(A_new[:nx*ny])
    {
        while (error > tolerance && itc < itc_max)
        {
            error = jacobi(A, A_new, nx, ny);
            swap(A, A_new, nx, ny);

            if(itc % 100 == 0) printf("%5d, %0.6f\n", itc, error);

            itc++;
        }
    }

    double end_time = omp_get_wtime();

    printf("Total run time = %f s.\n", end_time - start_time);

    free(A);
    free(A_new);

    return 0;
}


void initial(double *__restrict A, double *__restrict  A_new, int nx, int ny)
{
    memset(A, 0, nx * ny * sizeof(double));
    memset(A_new, 0, nx * ny * sizeof(double));

    for (int i = 0; i < nx; i++)
    {
        A[i] = 1.0;
        A_new[i] = 1.0;
    }
}

double jacobi(double *__restrict A, double *__restrict A_new, int nx, int ny)
{
    double error = 0.0;

    #pragma acc parallel loop reduction(max:error) copyin(A[:nx*ny]) copy(A_new[:nx*ny])
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            A_new[index(i, j, nx)] = 0.25 * (A[index(i-1, j, nx)] + A[index(i+1, j, nx)]
                                            + A[index(i, j-1, nx)] + A[index(i, j+1, nx)]);
            error = fmax(error, fabs(A_new[index(i, j, nx)] - A[index(i, j, nx)]));
        }
    }

    return error;   
}

void swap(double *__restrict A, double *__restrict A_new, int nx, int ny)
{
    #pragma acc parallel loop copyin(A_new[:nx*ny]) copyout(A[:nx*ny]) tile(1,1024)
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            A[index(i, j, nx)] = A_new[index(i, j, nx)];
        }
    }
}
