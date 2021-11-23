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
#include<mpi.h>

#define index(x, y, ny) ((x) * (ny) + (y))    // a transform of 2-d array index

double jacobi(double *A, double *A_new, int nx, int ny);

/*  data and block mesh arrangement
 *  ny
 *  ^
 *  | (0,2)  
 *  | (0,1)   
 *  | (0,0)  (1,0)  (2,0)
 *  ---------------------> nx
 * 
 *  for c and cpp language, data fill column first, (1,1),(1,2) ...
 *  rank assign : column first: rank0 = (1,1), rank2 = (1,2) ...
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

    /* number of blocks in x and y direction */
    int nx_block = 9;
    int ny_block = 1;

    /* wrong number of processer */
    if(nx_block * ny_block != num_process) return 1;

    /* block index */
    int block_x = (int) (rank / ny_block);
    int block_y = rank % ny_block;

    // divide by row
    int weight = total_nx / nx_block + 2;
    int height = total_ny / ny_block + 2;
    if (block_x == 0 || block_x == nx_block-1) // only have one boundary
    {
        weight = weight - 1;
    }
    if (block_y == 0 || block_y == ny_block-1) // only have one boundary
    {
        height = height - 1;
    }

    const int nx = weight;
    const int ny = height;
    
    double * A = (double *)malloc(sizeof(double) * nx * ny);
    double * A_new = (double *)malloc(sizeof(double) * nx * ny);

    double * sender_top = (double *)malloc(sizeof(double) * nx);
    double * sender_bottom = (double *)malloc(sizeof(double) * nx);
    // double * sender_left = (double *)malloc(sizeof(double) * ny);
    // double * sender_right = (double *)malloc(sizeof(double) * ny);

    double * receiver_top = (double *)malloc(sizeof(double) * nx);
    double * receiver_bottom = (double *)malloc(sizeof(double) * nx);
    // double * receiver_left = (double *)malloc(sizeof(double) * ny);
    // double * receiver_right = (double *)malloc(sizeof(double) * ny);

    // initial
    memset(A, 0, nx * ny * sizeof(double));
    memset(A_new, 0, nx * ny * sizeof(double));
    if(block_x == 0)
    {
        /* set left column = 1.0 */
        for (int i = 0; i < ny; i++)
        {
            A[i] = 1.0;
            A_new[i] = 1.0;
        }
    }

    if (rank == 0)
    {
        printf("Solve 2-D Laplace equation using jacobi relaxation\nmesh size: %d x %d\n", total_nx, total_ny);
        printf("Using %d cores.\n", num_process);
        printf("Block mesh: %d x %d\n", nx_block, ny_block);        
    }

    // printf("I am rank %d, block : %d x %d\n", rank, block_x, block_y);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // jacobi interation
    int itc = 0;
    double error, tmp;
    double error_max = 1.0;
    while (error_max > tolerance && itc < itc_max)
    {
        error = 0.0;
        MPI_Request req_tops, req_topr, req_bots, req_botr, req_lefs, req_lefr, req_rigs, req_rigr;

        /* calculate left column */
        for (int i = 1, j = 2; j < ny-2; j++)
        {
            A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                            + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
            tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
            error = error > tmp ? error : tmp;
        }

        /* exchange message with left
         * message tag: 0 - message from/to left; 1 - message from/to right (for sender) */
        if (block_x > 0)
        {
            MPI_Isend(A_new+ny+1, ny-2, MPI_DOUBLE, rank-ny_block, 0, MPI_COMM_WORLD, &req_lefs);
            MPI_Irecv(A_new+1, ny-2, MPI_DOUBLE, rank-ny_block, 1, MPI_COMM_WORLD, &req_lefr);
        }

        /* calculate right column */
        for (int i = nx-2, j = 2; j < ny-2; j++)
        {
            A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                            + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
            tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
            error = error > tmp ? error : tmp;
        }

        /* exchange message with right */
        if (block_x < nx_block-1)
        {
            MPI_Isend(A_new + (nx-2) * ny + 1, ny - 2, MPI_DOUBLE, rank + ny_block, 1, MPI_COMM_WORLD, &req_rigs);
            MPI_Irecv(A_new + (nx-1) * ny + 1, ny - 2, MPI_DOUBLE, rank + ny_block, 0, MPI_COMM_WORLD, &req_rigr);
        }

        /* calculate bottom row */
        for (int j = 1, i = 2; i < nx-2; i++)
        {
            A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                            + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
            tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
            error = error > tmp ? error : tmp;
        }

        /* exchange message with bottom */
        if (block_y > 0)
        {
            for (int i = 1; i < nx-1; i++)
            {
                sender_bottom[i] = A_new[index(i, 1, ny)];
            }
            MPI_Isend(sender_bottom+1, nx - 2, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_bots);
            MPI_Irecv(receiver_bottom+1, nx - 2, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &req_botr);
        }        

        /* calculate top row */
        for (int j = ny-2, i = 2; i < nx-2; i++)
        {
            A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                            + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
            tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
            error = error > tmp ? error : tmp;
        }

        /* exchange message with top
         * message tag: 0 - message from/to bottom; 1 - message from/to top (for sender) */
        if (block_y < ny_block-1)
        {
            for (int i = 1; i < nx-1; i++)
            {
                sender_top[i] = A_new[index(i, ny-2, ny)];
            }
            MPI_Isend(sender_top+1, nx - 2, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &req_tops);
            MPI_Irecv(receiver_top+1, nx - 2, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_topr);
        }

        /* calculate inner data */
        for (int i = 2; i < nx-2; i++)
        {
            for (int j = 2; j < ny-2; j++)
            {
                A_new[index(i, j, ny)] = 0.25 * (A[index(i-1, j, ny)] + A[index(i+1, j, ny)]
                                                + A[index(i, j-1, ny)] + A[index(i, j+1, ny)]);
                tmp = fabs(A_new[index(i, j, ny)] - A[index(i, j, ny)]);
                error = error > tmp ? error : tmp;
            }
        }

        /* ---------------- wait messages --------------- */
        /* wait messages of left */
        if (block_x > 0)
        {
            MPI_Wait(&req_lefs, MPI_STATUS_IGNORE);
            MPI_Wait(&req_lefr, MPI_STATUS_IGNORE);
        }

        /* wait messages of right */
        if (block_x < nx_block-1)
        {
            MPI_Wait(&req_rigs, MPI_STATUS_IGNORE);
            MPI_Wait(&req_rigr, MPI_STATUS_IGNORE);
        }

        /* wait messages of bottom */
        if (block_y > 0)
        {
            MPI_Wait(&req_bots, MPI_STATUS_IGNORE);
            MPI_Wait(&req_botr, MPI_STATUS_IGNORE);
            for (int i = 1; i < nx-1; i++)
            {
                A_new[index(i, 0, ny)] = receiver_bottom[i];
            }
        }

        /* wait messages of top */
        if (block_y < ny_block-1)
        {
            MPI_Wait(&req_tops, MPI_STATUS_IGNORE);
            MPI_Wait(&req_topr, MPI_STATUS_IGNORE);
            for (int i = 1; i < nx-1; i++)
            {
                A_new[index(i, ny-1, ny)] = receiver_top[i];
            }
        }


        for (int i = 1; i < nx-1; i++)
        {
            for (int j = 1; j < ny-1; j++)
            {
                A[index(i, j, ny)] = A_new[index(i, j, ny)];
            }
        }

        /* exchange message with corners  * message tag:     1   3
         *                                  (for senders)    0   2 */
        /* exchange message with top-right corner --- send then receive */
        if (block_x < nx_block-1 && block_y < ny_block-1)
        {
            MPI_Send(A + (nx-1)*ny - 2, 1, MPI_DOUBLE, rank + 1 + ny_block, 3, MPI_COMM_WORLD);
            MPI_Recv(A + nx*ny - 1, 1, MPI_DOUBLE, rank + 1 + ny_block, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* exchange message with bottom-left corner --- receive then send */
        if (block_x > 0 && block_y > 0)
        {
            MPI_Recv(A, 1, MPI_DOUBLE, rank - 1 - ny_block, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(A + ny + 1, 1, MPI_DOUBLE, rank - 1 - ny_block, 0, MPI_COMM_WORLD);
        }

        /* exchange message with top-left corner --- send then receive */
        if (block_x > 0 && block_y < ny_block-1)
        {
            MPI_Send(A + ny + ny - 2, 1, MPI_DOUBLE, rank - ny_block + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(A + ny - 1, 1, MPI_DOUBLE, rank - ny_block + 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* exchange message with bottom-right corner --- receive then send */
        if (block_x < nx_block-1 && block_y > 0)
        {
            MPI_Recv(A + (nx-1)*ny, 1, MPI_DOUBLE, rank + ny_block - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(A + (nx-2)*ny + 1, 1, MPI_DOUBLE, rank + ny_block - 1, 2, MPI_COMM_WORLD);
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