/* This code simulates a 2-dimensional laplace equation
 * using jacobi interation
 * accelated by mpi
 * divide data zone into 2d mesh
 * author: btli(2021)
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define index(x, y, ny) ((x) * (ny) + (y))    // a transform of 2-d array index

double jacobi(double *A, double *A_new, int nx, int ny);
void swap(double *A, double *A_new, int nx, int ny);

/*  data arrangement
 *  ny
 *  ^
 *  | (0,2)  
 *  | (0,1)   
 *  | (0,0)  (1,0)  (2,0)
 *  ---------------------> nx
 * 
 *  for c and cpp language, data fill column first, (0,0),(0,1) ...
*/
int main()
{
    const int nx = 4320;
    const int ny = 4320;
    const int itc_max = 1000;

    const double tolerance = 1e-5;

    
    double * A = (double *)malloc(sizeof(double) * nx * ny);
    double * A_new = (double *)malloc(sizeof(double) * nx * ny);


    // initial
    memset(A, 0, nx * ny * sizeof(double));
    memset(A_new, 0, nx * ny * sizeof(double));

    /* set top boundary = 1.0 */
    for (int i = 0; i < nx; i++)
    {
        A[index(i, ny-1, ny)] = 1.0;
        A_new[index(i, ny-1, ny)] = 1.0;
    }

    // jacobi interation
    int itc = 0;
    double error = 1.0;
    while (error > tolerance && itc++ < itc_max)
    {
        error = jacobi(A, A_new, nx, ny);
        swap(A, A_new, nx, ny);

        if(itc % 100 == 0) printf("%5d, %0.6f\n", itc, error);
    }    

    free(A);
    free(A_new);

    return 0;
}



double jacobi(double *A, double *A_new, int nx, int ny)
{
    double tmp, error = 0.0;
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                            + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
            tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
            error = error > tmp ? error : tmp;
        }
    }

    return error;
}

void swap(double *A, double *A_new, int nx, int ny)
{
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            A[index(i, j, ny)] = A_new[index(i, j, ny)];
        }
    }
}