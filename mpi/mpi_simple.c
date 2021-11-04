#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char **argv)
{
    int size = 4096;

    double start, stop;
    int i, j, k, m, error = 0;
    double *a, *b, *c;
    int rank, numprocs, line;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("We are using the processor %s, rank %d from total %d processor.\n", processor_name, rank, numprocs);

    line = size / numprocs;
    
    // main process
    if (rank == 0)
    {     
        // initialize the random matrix
        a = (double *)malloc(sizeof(double) * size * size);
        b = (double *)malloc(sizeof(double) * size * size);
        c = (double *)malloc(sizeof(double) * size * size);
        srand((unsigned)time(NULL));
        for (i = 0; i < size * size; i++)
        {
            a[i] = (float)rand() / (RAND_MAX);
            b[i] = (float)rand() / (RAND_MAX);
            c[i] = 0.0;
        }

        start = MPI_Wtime(); // start conting

        // send matrix b to each process
        for (i = 1; i < numprocs; i++)
        {
            MPI_Send(b, size * size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        // send the row block of matrix a to each process
        for (m = 1; m < numprocs; m++)
        {
            MPI_Send(a + m * line * size, size * line, MPI_DOUBLE, m, 1, MPI_COMM_WORLD);
        }

        // main process compute the first row block
        for (i = 0; i < line; i++)  // only comput first row block
        {
            for (j = 0; j < size; j++)
            {           
                for (k = 0; k < size; k++)
                {
                    c[i * size + j] += a[i * size + k] * b[k * size + j];
                }
            }
        }

        // receive the results from other process
        for (k = 1; k < numprocs; k++)
        {
            MPI_Recv(c + k * line * size, line * size, MPI_DOUBLE, k, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // count time
        stop = MPI_Wtime();

        printf("rank:%d time:%lfs\n", size, stop - start);

        free(a);
        free(b);
        free(c);

    }
    else    //other process
    {
        a = (double *)malloc(sizeof(double) * line * size);   // row block of matrix a
        b = (double *)malloc(sizeof(double) * size * size);
        c = (double *)malloc(sizeof(double) * line * size);   // row block of matrix c
        // receive matrix b
        MPI_Recv(b, size * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // receive the row block of matrix a
        MPI_Recv(a, size * line, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //compute the row of matrix c
        for (i = 0; i < line; i++)// only comput first row block
        {
            for (j = 0; j < size; j++)
            {
                c[i * size + j] = 0;
                for (k = 0; k < size; k++)
                {
                    c[i * size + j] += a[i * size + k] * b[k * size + j];
                }
            }
        }
        // send the result to main process
        MPI_Send(c, line * size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        free(a);
        free(b);
        free(c);
    }

    MPI_Finalize(); 

    return 0;
}