#include <stdio.h>
#include <omp.h>
void pooh(int id, double * A);


int main()
{
    double A[30] = {0};
    int idd, size;
    double *B;

    omp_set_num_threads(12);

    #pragma omp parallel
    {
        double *C;
        int id = omp_get_thread_num();
        int Ntds = omp_get_num_threads();
        if(id == 0) size = Ntds;
        pooh(id, A);
        
        printf("A of ID(%d) = %lf \n", id, A[id]);
        if(id % 2 == 0) idd = id;
        printf("idd = %d, \n",idd);
        // idd = idd + 1;
        // printf("idd = %d\n", idd);

        printf("B = %x, B = %x \n", B, C);
    }
    printf("%d threads.", size);
    return 0;
}

void pooh(int id, double * A)
{
    A[id] = id;
}