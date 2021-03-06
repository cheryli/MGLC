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

#define index(x, y, ny) ((x) * (ny) + (y))    // a transform of 2-d array index

double jacobi(double *A, double *A_new, int nx, int ny);
void swap(double *A, double *A_new, int nx, int ny);

/*  data mesh arrangement
 *  ny
 *  ^
 *  | (0,2)  
 *  | (0,1)   
 *  | (0,0)  (1,0)  (2,0)
 *  ---------------------> nx
 * 
    for c and cpp language, data fill column first, (1,1),(1,2) ...
*/


int main()
{
    const int total_nx = 4320;
    const int total_ny = 4320;
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

    // divide by column
    int weight = total_nx / num_process + 2;
    if (rank == 0 || rank == num_process-1) // only have one boundary for message exchange
    {
        weight--;
    }

    // handle 'Mesh not divisible'
    if (rank < total_nx % num_process)
    {
        weight++;
    }

    const int nx = weight;
    const int ny = total_ny;
    
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

    if(rank == 0)
    {
        printf("Solve 2-D Laplace equation using jacobi relaxation\nmesh size: %d x %d\n", total_nx, total_ny);
        printf("Using %d cores.", num_process);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // jacobi interation
    int itc = 0;
    double error;
    double error_max = 1.0;
    while (error_max > tolerance && itc++ < itc_max)
    {
        error = jacobi(A, A_new, nx, ny);
        swap(A, A_new, nx, ny);

        if(rank > 0)
        {
            // exchange message with top --- receive then send
            MPI_Recv(A, ny, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(A + ny, ny, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank < num_process-1)
        {
            // exchange message with bottom --- send then receive
            MPI_Send(A + (nx-2)*ny, ny, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
            MPI_Recv(A + (nx-1)*ny, ny, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&error, &error_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if(rank == 0)
        {
            if(itc % 100 == 0) printf("%5d, %0.6f\n", itc, error_max);
        }
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