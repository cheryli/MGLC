! /* This code simulates a 2-dimensional laplace equation
!  * using jacobi interation
!  * accelated by mpi
!  * author: btli(2021)
!  */

! /*  data mesh arrangement
!  *  ny
!  *  ^
!  *  | (1,3)  
!  *  | (1,2)   
!  *  | (1,1)  (2,1)  (3,1)
!  *  ---------------------> nx
!  * 
!     for fortran language, data fill row first, (1,1),(2,1) ...
! */

program main
    use mpi_f08
    use jacobi
    implicit none
    
    integer, parameter :: fp_kind=kind(1.0d0)
    integer, parameter :: total_nx = 4500, total_ny = 4500, itc_max = 1000
    real(fp_kind), parameter :: tolerance = 1.0e-5_fp_kind
    integer :: nx, ny
    integer :: i, j, itc
    real(fp_kind),allocatable   :: A(:, :)
    real(fp_kind),allocatable   :: A_new(:, :)
    real(fp_kind) :: error, error_max
    real(fp_kind) :: start_time, end_time
    integer :: rc, rank, num_process, name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    ! type(mpi_status) :: status

    call mpi_init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call mpi_get_processor_name(processor_name, name_len, rc)

    ! divide by row
    nx = total_nx
    ny = total_ny / num_process + 2;
    if (rank == 0 .or. rank == num_process-1) then !! only have one boundary
        ny = ny - 1
    endif

    allocate ( A(0:nx-1,0:ny-1), A_new(0:nx-1,0:ny-1) )

    !! initial
    A = 0.0d0
    A_new = 0.0d0

    
    if (rank == 0) then
        A(:, 0) = 1.0d0
        A_new(:, 0) = 1.0d0
    endif


    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()


    !! jacobi iteration
    itc = 0
    error_max = 1.0
    do while (error_max > tolerance .and. itc < itc_max)
        error = calc_next(A, A_new, nx, ny)
        call swap(A, A_new, nx, ny)

        if(rank > 0) then
            call MPI_Recv(A, nx, MPI_DOUBLE_PRECISION, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
            call MPI_Send(A(0,1), nx, MPI_DOUBLE_PRECISION, rank-1, 0, MPI_COMM_WORLD, rc)
        endif

        if(rank < num_process-1) then
            call MPI_Send(A(0, ny-2), nx, MPI_DOUBLE_PRECISION, rank+1, 1, MPI_COMM_WORLD, rc)
            call MPI_Recv(A(0, ny-1), nx, MPI_DOUBLE_PRECISION, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
            
        endif

        call MPI_Barrier(MPI_COMM_WORLD, rc)
        call MPI_Allreduce(error, error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
        if(rank == 0) then
            if(MOD(itc, 100) == 0) then
                write(*, *) itc, error_max
            endif
        endif

        itc = itc + 1
    enddo


    call MPI_Barrier(MPI_COMM_WORLD, rc);
    end_time = MPI_Wtime()
    if (rank == 0) then
        write(*,*) "Total run time =", end_time-start_time, "s."
    endif

    deallocate (A, A_new)
    call MPI_Finalize(rc)

end program main


