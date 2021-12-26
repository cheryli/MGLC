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
    use mpi
    implicit none
    
    integer, parameter :: total_nx = 2000, total_ny = 2000, itc_max = 1000
    real(8), parameter :: tolerance = 1.0e-5
    integer :: nx, ny
    integer :: i, j, itc
    real(8),allocatable   :: A(:, :), A_new(:, :), f(:, :), A_p(:, :)
    real(8) :: error, error_max, check_diff
    real(8) :: start_time, end_time
    integer :: rc, rank, num_process, name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
    integer :: dims(2) = (0, 0), coords(2)
    logical :: periods(2)
    data periods/2*.false./
    integer :: comm2d, rank2d, column_y
    integer :: nbr_left, nbr_right, nbr_top, nbr_bottom
    

    ! type(mpi_status) :: status

    call mpi_init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)


    !!! ---- Decomposition the domain
    ! dims(1) = 9
    ! dims(2) = 1
    call MPI_Dims_create(num_process, 2, dims, rc)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, rc)
    if(rank == 0) then
        write(*,*) "dimens is x*y = ", dims(1), "x", dims(2)
    endif

    ! get rank in decomposition
    call MPI_Comm_rank(comm2d, rank2d, rc)
    ! write(*,*) "process ", rank, " of total ", num_process, "is alive."

    ! get the neighbors
    call MPI_Cart_shift(comm2d, 0, 1, nbr_left, nbr_right, rc)
    call MPI_Cart_shift(comm2d, 1, 1, nbr_bottom, nbr_top, rc)
    ! write(*,*) "I'm process ", rank, "My neighbor lrbt is", nbr_left, nbr_right, nbr_bottom, nbr_top
    
    ! determain sub-domain size
    call MPI_Cart_get(comm2d, 2, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(1), dims(1))
    call decompose_1d(total_ny, ny, coords(2), dims(2))
    ! if(rank == 0) then
    !     write(*,*) "dims = ", dims(1), dims(2)
    ! endif   
    ! write(*,*) "coords = ", coords(1), coords(2)
    ! write(*,*) "nx*ny = ", nx, ny

    ! allocate a ghost layer at each boundary
    allocate(A(0:nx+1, 0:ny+1))
    allocate(A_new(0:nx+1, 0:ny+1))
    allocate(f(0:nx+1, 0:ny+1)) 
    allocate(A_p(0:nx+1, 0:ny+1))

    ! construct the datatype for the exchange in y direction(non-contiguous)
    call MPI_Type_vector(ny+2, 1, nx+2, MPI_DOUBLE_PRECISION, column_y, rc)     ! two ghost layers
    call MPI_Type_commit(column_y, rc)

    call init(A, A_new, f, nx, ny, dims, coords)

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    !! jacobi iteration
    itc = 0
    error_max = 1.0
    A_p = A
    do while(error_max > tolerance .AND. itc <= itc_max)
        itc =  itc + 2
        call jacobi_with_message_exchange(A, A_new, f, nx, ny, coords, comm2d, column_y, &
                                nbr_left, nbr_right, nbr_bottom, nbr_top)

        call jacobi_with_message_exchange(A_new, A, f, nx, ny, coords, comm2d, column_y, &
                                nbr_left, nbr_right, nbr_bottom, nbr_top)

        if(MOD(itc, 100) == 0) then
            error = check_diff(A, A_p, nx, ny)
            call MPI_Allreduce(error, error_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm2d, rc)
            if(rank == 0) then
                write(*, *) itc, error_max
            endif
        endif
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc);
    end_time = MPI_Wtime()
    if (rank == 0) then
        write(*,*) "Total run time =", end_time-start_time, "s."
    endif
    
    deallocate(A)
    deallocate(A_new)
    deallocate(A_p)
    deallocate(f)

    call MPI_Finalize(rc)

end program main


subroutine decompose_1d(total_n, local_n, rank, num_process)
    implicit none
    integer, intent(in) :: total_n, rank, num_process
    integer, intent(out) :: local_n

    local_n = total_n / num_process

    if (rank < MOD(total_n, num_process)) then
        local_n = local_n + 1
    endif

end subroutine decompose_1d


subroutine init(A, A_new, f, nx, ny, dims, coords)
    implicit none
    integer, intent(in) :: nx, ny, dims(2), coords(2)
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1), f(0:nx+1, 0:ny+1)
    integer :: i, j

    ! initialize
    do j = 0, ny+1
        do i = 0, nx+1
            A(i, j) = 0.0d0
            A_new(i, j) = 0.0d0
            f(i, j) = 0.0d0
        enddo
    enddo

    ! boundary condition (top boundary = 1.0, others = 0.0)
    if (coords(2) == dims(2) - 1) then
        do i = 0, nx+1 
            A(i, ny+1) = 1.0d0
            A_new(i, ny+1) = 1.0d0
        enddo
    endif

end subroutine init


subroutine jacobi_with_message_exchange(A, A_new, f, nx, ny, coord, comm2d, column_y, nbr_left, nbr_right, nbr_bottom, nbr_top)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, comm2d, column_y, coord(2)
    integer, intent(in) :: nbr_left, nbr_right, nbr_top, nbr_bottom
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1), f(0:nx+1, 0:ny+1)
    integer :: i, j, rc, req(8)

    ! start message change at boundary
    ! corners don't matter in this problem

    ! exchange messages along y (ghost row)
    ! message passing to top(j++)
    call MPI_Irecv(A(1, 0), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 1, comm2d, req(1), rc)
    call MPI_Isend(A(1, ny), nx, MPI_DOUBLE_PRECISION, nbr_top, 1, comm2d, req(2), rc)
    
    ! message passing to bottom(j--)
    call MPI_Irecv(A(1, ny+1), nx, MPI_DOUBLE_PRECISION, nbr_top, 0, comm2d, req(3), rc)
    call MPI_Isend(A(1, 1), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 0, comm2d, req(4), rc)

    ! exchange messages along x (ghost column)
    ! message passing to right(i++)
    call MPI_Irecv(A(0, 1), 1, column_y, nbr_left, 0, comm2d, req(5), rc)
    call MPI_Isend(A(nx, 1), 1, column_y, nbr_right, 0, comm2d, req(6), rc)
    
    ! message passing to left(i--)
    call MPI_Irecv(A(nx+1, 1), 1, column_y, nbr_right, 1, comm2d, req(7), rc)
    call MPI_Isend(A(1, 1), 1, column_y, nbr_left, 1, comm2d, req(8), rc)

    ! update inner points
    do j = 2, ny-1
        do i = 2, nx-1
            A_new(i, j) = 0.25 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) + f(i, j)) 
        enddo
    enddo

    ! wait for message exchange finish
    call MPI_Waitall(8, req, MPI_STATUSES_IGNORE, rc)

    ! update the boundary
    do j = 1, ny, ny-1
        do i = 1, nx
            A_new(i, j) = 0.25 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) + f(i, j)) 
        enddo
    enddo
    do i = 1, nx, nx-1
        do j = 1, ny
            A_new(i, j) = 0.25 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) + f(i, j)) 
        enddo
    enddo

end subroutine jacobi_with_message_exchange


subroutine exchange_message_nonblock(A, nx, ny, coord, comm2d, column_y, nbr_left, nbr_right, nbr_bottom, nbr_top)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, comm2d, column_y, coord(2)
    integer, intent(in) :: nbr_left, nbr_right, nbr_top, nbr_bottom
    real(8), intent(in) :: A(0:nx+1, 0:ny+1)
    integer :: i, j, rc, req(8)

    ! corners don't matter in this problem

    ! exchange messages along y (ghost row)
    ! message passing to top(j++)
    call MPI_Irecv(A(1, 0), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 1, comm2d, req(1), rc)
    call MPI_Isend(A(1, ny), nx, MPI_DOUBLE_PRECISION, nbr_top, 1, comm2d, req(2), rc)
    
    ! message passing to bottom(j--)
    call MPI_Irecv(A(1, ny+1), nx, MPI_DOUBLE_PRECISION, nbr_top, 0, comm2d, req(3), rc)
    call MPI_Isend(A(1, 1), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 0, comm2d, req(4), rc)

    ! exchange messages along x (ghost column)
    ! message passing to right(i++)
    call MPI_Irecv(A(0, 1), 1, column_y, nbr_left, 0, comm2d, req(5), rc)
    call MPI_Isend(A(nx, 1), 1, column_y, nbr_right, 0, comm2d, req(6), rc)
    
    ! message passing to left(i--)
    call MPI_Irecv(A(nx+1, 1), 1, column_y, nbr_right, 1, comm2d, req(6), rc)
    call MPI_Isend(A(1, 1), 1, column_y, nbr_left, 1, comm2d, req(7), rc)

    call MPI_Waitall(8, req, MPI_STATUSES_IGNORE, rc)

end subroutine exchange_message_nonblock





function check_diff(A, A_p, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: A(0:nx+1, 0:ny+1)
    real(8), intent(inout) :: A_p(0:nx+1, 0:ny+1)
    integer :: i, j
    real(8) :: error, check_diff

    error = 0.0d0
    do j = 1, ny
        do i = 1, nx
            error = max( error, abs(A_p(i,j)-A(i,j)) )
        enddo
    enddo

    A_p = A
    
    check_diff = error

end function check_diff


subroutine swap(A, A_new, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1)
    integer :: i, j

    do j = 1, ny
        do i = 1, nx
            A(i, j) = A_new(i, j)
        enddo
    enddo

end subroutine swap