---
layout: default
title: "MPI parallel of Jacobi iteration"
permalink: /mpi/jacobi/
---


## Jacobi iteration of Laplace equation.
We can use jacobi iteration to numerically solve the Laplace's equation.  
The general form of Laplace's equation is:

$$
\begin{equation}
    \nabla^2 u = 0
\end{equation}
$$

For simplicity, we consider the two-dimensional case in Cartesian coordinates, it takes the form:

$$
\begin{equation}
    \left(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}\right) u(x, y) = 0
\end{equation}
$$

And then we use the uniform mesh and replace the partial differences with the second-order center difference scheme on the mesh as following:

$$
\begin{equation}
    \left(\frac{\partial^2 u}{\partial x^2}\right)_{i, j} = \frac{u_{i+1, j} - 2u_{i, j} + u_{i-1, j}}{(\Delta x)^2} 
\end{equation}
$$

$$
\begin{equation}
    \left(\frac{\partial^2 u}{\partial y^2}\right)_{i, j} = \frac{u_{i, j+1} - 2u_{i, j} + u_{i, j-1}}{(\Delta y)^2} 
\end{equation}
$$

Suppose $\Delta x = \Delta y = 1$, then we can rewrite the Poisson's equation as:

$$
\begin{equation}
    4u_{i, j} - u_{i-1, j} - u_{i+1, j} - u_{i, j-1} - u_{i, j+1} = 0
\end{equation}
$$

Remember $i, j$ could be any point in the computation domain. So above equation is actually a large sparse systems of equations, and can be solved iteratively using jacobi method:

$$
\begin{equation}
    u_{i,j}^{n+1} = \left(\frac{u_{i-1, j} + u_{i+1, j} + u_{i, j-1} + u_{i, j+1}}{4} \right)^n
\end{equation}
$$

where $n$ is the iteration step. This equation is simple and clear, at every new step, the value of a mesh point is the average of its surrounding points at last step.

In study of heat conduction($u = T$), Laplace equation describe a steady state temperature distribution. If we set the temperature at upper boundary equals 1, and the temperature at other boundary equals 0, the temperature distribution at steady state looks like this:

![laplace-steady-state](/assets/laplace.jpg)

## Code structure of serial program
It's quite a simple problem, you can view the serial code [here](). 

We use two arrays namely `A` and `A_new` to store the temperature matrix. And for simplicity, we define a macro in C language to index a two dimension array.
```c
#define index(x, y, ny) ((x) * (ny) + (y))    // a transform of 2-d array index
```
Then the main iteration process is putted into a single `jacobi` function, and after each iteration step, we need to swap `A` and `A_new` by `swap` function:
```c
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
```

And what the `main` function do is simply `initial` the value, and then do the iteration until reach to its convergence or maximum iteration step.
```c
while (error_max > tolerance && itc++ < itc_max)
{
    error = jacobi(A, A_new, nx, ny);
    swap(A, A_new, nx, ny);

    if(itc % 100 == 0) printf("%5d, %0.6f\n", itc, error);
}
```

## MPI parallelization of jacobi iteration
One way to accelerate above serial program on multi-processor systems is MPI(Message Passing Interface). And we will introduce the block division and latency hiding of this MPI implementation.

we will start with the blocked version(without latency hiding). The idea is, as illustration,we can divide the computation domain into many sub-domains, and each processor will handle one sub-domain. The update of inner points is easy and totally local, but the update of boundary points acquires the information of other sub-domains. So we need a layer of ghost points to exchange messages after every steps of iteration.

![block_division_1d](/assets/block_division_1d.jpg)

But you can see, as we increase the number of processor, the computation cost of each processor will reduce but the message passing cost remain unchanged. It's certainly not a good choice for massive parallel tasks. A better way is use the 2d division:

![block_division_1d](/assets/block_division_2d.jpg)


###  Block index
This 2d division method reduce the message passing cost, but increase the program complexity. We give each sub-domain a block index, and assign the rank along y first.
```c
/*  data and block mesh arrangement
 *  ny
 *  ^
 *  | (0,2)  
 *  | (0,1)   
 *  | (0,0)  (1,0)  (2,0)
 *  ---------------------> nx
 * 
 *  for c and cpp language, data fill column first, (0,0),(0,1) ...
 *  rank assign : column first: rank0 = (0,0), rank2 = (1,0) ...
*/
```

Before the iteration, we calculate the block index first:
```c
/* block index */
int block_x = (int) (rank / ny_block);
int block_y = rank % ny_block;
```
In order to set the top boundary of origin domain, we now set the top boundary of domain those `block_y = ny_block` equal 1.
```c
if(block_y == ny_block)
{
    /* set top boundary = 1.0 */
    for (int i = 0; i < ny; i++)
    {
        A[index(i, ny-1, ny)] = 1.0;
        A_new[index(i, ny-1, ny)] = 1.0;
    }
}
```

### Message exchange
And the message passing becomes more complicated. For each sub-domain, for example, we need to exchange message with its top boundary. We now use its "block index" to identify if it has the top sub-domain and the rank id of its top sub-domain. And because the points at top boundary are discontinuous in virtual memory, we need a buffer `sender_tb` and `receiver_tb`:
```c
/* exchange message with top -- send then receive
  * message tag: 0 - message from/to bottom; 1 - message from/to top (for sender) */
if (block_y < ny_block-1)
{
    for (int i = 1; i < nx-1; i++)
    {
        sender_tb[i] = A[index(i, ny-2, ny)];
    }
    MPI_Send(sender_tb+1, nx - 2, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
    MPI_Recv(receiver_tb+1, nx - 2, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 1; i < nx-1; i++)
    {
        A[index(i, ny-1, ny)] = receiver_tb[i];
    }
}
```
Recall the sketch map of block 2d-division. we don't need to exchange the corner message with its top side sub-domain. That's the biggest difference between 1-d and 2-d division. We need to handle corner message especially. For example, we exchange message with top-right corner as follows.
```c
/* exchange message with top-right corner --- send then receive */
if (block_x < nx_block-1 && block_y < ny_block-1)
{
    MPI_Send(A + (nx-1)*ny - 2, 1, MPI_DOUBLE, rank + 1 + ny_block, 3, MPI_COMM_WORLD);
    MPI_Recv(A + nx*ny - 1, 1, MPI_DOUBLE, rank + 1 + ny_block, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
```

### Mesh not divisible ?
If `nx = 100` and we want to set `nx_block = 3`, what should we do? The answer is simple. We can set first column block `nx = 34`, and all other block `nx = 33`, by following codes:
```c
// divide by column
int weight = total_nx / nx_block + 2;
int height = total_ny / ny_block + 2;
if (block_x == 0 || block_x == nx_block-1) // only have one boundary
{
    weight--;
}
if (block_y == 0 || block_y == ny_block-1) // only have one boundary
{
    height--;
}

// handle 'Mesh not divisible'
if (block_x < total_nx % nx_block)
{
    weight++;
}
if (block_y < total_ny % ny_block)
{
    height++;
}

const int nx = weight;
const int ny = height;
```

### Latency hiding by non-block message passing

The latency of message passing is the worst
enemy of parallel efficiency when using MPI. A important improvement we can make is to overlap the computation and message passing to hide the latency. So we need to use the non-block message passing function.

For this simple problem, we can calculate the boundary points first, and then start sending the messages using non-block sending functions, take the left column as an example:
```c
// calculate left column first
error = jacobi(A, A_new, 3, ny);

// exchange message with left
if (rank > 0)
{
    MPI_Isend(A_new + ny, ny, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqts);
    MPI_Irecv(A_new, ny, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqtr);  
}
```
After all non-block sending start, we then  start calculate the inner points:
```c
// calculate middle columns finally
tmp = jacobi(A + ny, A_new + ny, nx-1, ny);
error = error > tmp ? error : tmp;
```
Before moving to the next step, we need to make sure the message passing has finished:
```c
// wait message of left
if (rank > 0)
{
    MPI_Wait(&reqts, MPI_STATUS_IGNORE);
    MPI_Wait(&reqtr, MPI_STATUS_IGNORE);
}
```






test


![this is a link to a wallpaper](/assets/wallpaper-1.jpg)


