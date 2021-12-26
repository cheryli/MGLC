subroutine exchange_message(A, nx, ny, coords, comm2d, column_y, nbr_left, nbr_right, nbr_bottom, nbr_top)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, comm2d, column_y, coords
    integer, intent(in) :: nbr_left, nbr_right, nbr_top, nbr_bottom
    real(8), intent(in) :: A(0:nx+1, 0:ny+1)
    integer :: i, j, rc

    ! corners don't matter in this problem

    ! exchange messages along y (ghost row)
    ! message passing to top(j++)
    call MPI_Sendrecv(A(1, ny), nx, MPI_DOUBLE_PRECISION, nbr_top, 0, &
                    A(1, 0), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 0, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    ! message passing to bottom(j--)
    call MPI_Sendrecv(A(1, 1), nx, MPI_DOUBLE_PRECISION, nbr_bottom, 1, &
                    A(1, ny+1), nx, MPI_DOUBLE_PRECISION, nbr_top, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange messages along x (ghost column)
    ! message passing to right(i++)
    call MPI_Sendrecv(A(nx, 1), 1, column_y, nbr_right, 0, &
                    A(0, 1), 1, column_y, nbr_left, 0, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    ! message passing to left(i--)
    call MPI_Sendrecv(A(1, 1), 1, column_y, nbr_left, 1, &
                    A(nx+1, 1), 1, column_y, nbr_right, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)


end subroutine exchange_message
