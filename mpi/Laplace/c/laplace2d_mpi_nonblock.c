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
    const int total_nx = 4500;
    const int total_ny = 4500;
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
    int weight = total_nx / num_process + 2;
    if (rank == 0 || rank == num_process-1) // only have one boundary
    {
        weight--;
    }

    const int nx = weight;
    const int ny = total_ny;
    
    double * A = (double *)malloc(sizeof(double) * nx * ny);
    double * A_new = (double *)malloc(sizeof(double) * nx * ny);

    // initial
    memset(A, 0, nx * ny * sizeof(double));
    memset(A_new, 0, nx * ny * sizeof(double));
    if(rank == 0)
    {
        for (int i = 0; i < ny; i++)
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
        error = 0.0;
        MPI_Request reqts, reqtr, reqbs, reqbr;
        if (rank == 0)
        {
            // receive message from right
            MPI_Irecv(A_new + (nx-1) * ny, ny, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &reqbr);
            for (int i = nx-2, j = 1; j < ny-1; j++)
            {
                A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                error = error > tmp ? error : tmp;
            }

            // send message to right
            // message tag 0 - message from/to top ; 1 - message from/to bottom 
            // message tag top/bottom for sender
            MPI_Isend(A_new + (nx-2) * ny, ny, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &reqbs);


            for (int i = 1; i < nx-2; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                    + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                    tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                    error = error > tmp ? error : tmp;
                }
            }

            MPI_Wait(&reqbs, MPI_STATUS_IGNORE);
            MPI_Wait(&reqbr, MPI_STATUS_IGNORE);          
        }
        else if (rank == num_process - 1)
        {
            // receive message from left
            MPI_Irecv(A_new, ny, MPI_DOUBLE, num_process-2, 1, MPI_COMM_WORLD, &reqtr);
            for (int i = 1, j = 1; j < ny-1; j++)
            {
                A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                error = error > tmp ? error : tmp;
            }

            // send message to left
            MPI_Isend(A_new + ny, ny, MPI_DOUBLE, num_process-2, 0, MPI_COMM_WORLD, &reqts);
            

            for (int i = 2; i < nx-1; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                    + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                    tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                    error = error > tmp ? error : tmp;
                }
            }

            MPI_Wait(&reqts, MPI_STATUS_IGNORE);
            MPI_Wait(&reqtr, MPI_STATUS_IGNORE);
        }
        else
        {
            /* receive message first */
            MPI_Irecv(A_new, ny, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqtr);
            MPI_Irecv(A_new + (nx-1)*ny, ny, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &reqbr);

            /* left column */
            for (int i = 1, j = 1; j < ny-1; j++)
            {
                A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                error = error > tmp ? error : tmp;
            }

            /* send message to left */
            MPI_Isend(A_new + ny, ny, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqts);
            

            /* right column */
            for (int i = nx-2, j = 1; j < ny-1; j++)
            {
                A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                error = error > tmp ? error : tmp;
            }

            /* send message to right */
            MPI_Isend(A_new + (nx-2)*ny, ny, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &reqbs);


            /* middle columns */
            for (int i = 2; i < nx-2; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                    + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                    tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                    error = error > tmp ? error : tmp;
                }
            }

            MPI_Wait(&reqts, MPI_STATUS_IGNORE);
            MPI_Wait(&reqtr, MPI_STATUS_IGNORE);
            MPI_Wait(&reqbs, MPI_STATUS_IGNORE);
            MPI_Wait(&reqbr, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i < nx-1; i++)
        {
            for (int j = 1; j < ny-1; j++)
            {
                A[index(i, j, ny)] = A_new[index(i, j, ny)];
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
