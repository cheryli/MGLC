subroutine exchange_message_nonblock(A, nx, ny, coords, comm2d, column_y, nbr_left, nbr_right, nbr_bottom, nbr_top)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, comm2d, column_y, coords(2)
    integer, intent(in) :: nbr_left, nbr_right, nbr_top, nbr_bottom
    real(8), intent(in) :: A(0:nx+1, 0:ny+1)
    integer :: rc, req(8)

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

    call MPI_Waitall(8, req, MPI_STATUSES_IGNORE, rc)

end subroutine exchange_message_nonblock
