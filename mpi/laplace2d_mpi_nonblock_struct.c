/* This code simulates a 2-dimensional laplace equation
 * using jacobi interation
 * accelated by mpi
 * author: btli(2021)
 */

#include<stdio.h> 
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

#define index(x, y, nx) ((x) * (nx) + (y))    // a transform of 2-d array index

double jacobi(double *A, double *A_new, int nx, int ny);

int main()
{
    const int total_nx = 4096;
    const int total_ny = 4096;
    const int itc_max = 1000;

    const double tolerance = 1e-5;


    MPI_Init(NULL, NULL);
    int num_process;
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // divide by row
    int height = total_ny / num_process + 2;
    if (rank == 0 || rank == num_process-1) // only have one boundary
    {
        height = total_ny / num_process + 1;
    }

    const int nx = total_nx;
    const int ny = height;
    
    double * A = (double *)malloc(sizeof(double) * nx * ny);
    double * A_new = (double *)malloc(sizeof(double) * nx * ny);

    // initial
    memset(A, 0, nx * ny * sizeof(double));
    memset(A_new, 0, nx * ny * sizeof(double));
    if(rank == 0)
    {
        for (int i = 0; i < nx; i++)
        {
            A[i] = 1.0;
            A_new[i] = 1.0;
        }
        printf("Solve 2-D Laplace equation using jacobi relaxation\nmesh size: %d x %d\n", total_nx, total_ny);
        printf("Using %d cores.", num_process);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // jacobi interation
    int itc = 0;
    double error, tmp;
    double error_max = 1.0;
    while (error_max > tolerance && itc < itc_max)
    {
        MPI_Request reqts, reqtr, reqbs, reqbr;

        // calculate top row first
        error = jacobi(A, A_new, nx, 3);

        // exchange message with top
        if (rank > 0)
        {
            MPI_Isend(A_new + nx, nx, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqts);
            MPI_Irecv(A_new, nx, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqtr);  
        }

        // calculate bottom row then
        tmp = jacobi(A + nx*(ny-3), A_new + nx*(ny-3), nx, 3);
        error = error > tmp ? error : tmp;

        /* exchange message with bottom */
        if (rank < num_process-1)
        {
            MPI_Isend(A_new + (ny-2)*nx, nx, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &reqbs);
            MPI_Irecv(A_new + (ny-1)*nx, nx, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &reqbr);  
        }

        // calculate middle rows finally
        tmp = jacobi(A + nx, A_new + nx, nx, ny-1);
        error = error > tmp ? error : tmp;

        
        // wait message of top
        if (rank > 0)
        {
            MPI_Wait(&reqts, MPI_STATUS_IGNORE);
            MPI_Wait(&reqtr, MPI_STATUS_IGNORE);
        }
        // wait message of bottom
        if (rank < num_process-1)
        {
            MPI_Wait(&reqbs, MPI_STATUS_IGNORE);
            MPI_Wait(&reqbr, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i < ny-1; i++)
        {
            for (int j = 1; j < nx-1; j++)
            {
                A[index(i, j, nx)] = A_new[index(i, j, nx)];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&error, &error_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if(rank == 0)
        {
            if(itc % 100 == 0) printf("%5d, %0.6f\n", itc, error_max);
        }
            
        itc++;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if(rank == 0) printf("Total run time = %f s.\n", end_time - start_time);
    

    free(A);
    free(A_new);

    MPI_Finalize(); 

    return 0;
}


double jacobi(double *A, double *A_new, int nx, int ny)
{
    double tmp, error = 0.0;
    for (int i = 1; i < ny-1; i++)
    {
        for (int j = 1; j < nx-1; j++)
        {
            A_new[index(i, j, nx)] = 0.25 * (A[index(i-1, j, nx)] + A[index(i+1, j, nx)]
                                            + A[index(i, j-1, nx)] + A[index(i, j+1, nx)]);
            tmp = fabs(A_new[index(i, j, nx)] - A[index(i, j, nx)]);
            error = error > tmp ? error : tmp;
        }
    }

    return error;
}