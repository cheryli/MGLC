subroutine message_passing_sendrecv()
    use mpi
    use commondata
    implicit none

    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Sendrecv(f_post(2, 1, ny), 1, row_x, nbr_top, 2, &
                    f_post(2, 1, 0), 1, row_x, nbr_bottom, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, 1, ny), 1, row_x, nbr_top, 5, &
                    f_post(5, 1, 0), 1, row_x, nbr_bottom, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(6, 1, ny), 1, row_x, nbr_top, 6, &
                    f_post(6, 1, 0), 1, row_x, nbr_bottom, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to bottom(j--)
    call MPI_Sendrecv(f_post(4, 1, 1), 1, row_x, nbr_bottom, 4, &
                    f_post(4, 1, ny+1), 1, row_x, nbr_top, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(7, 1, 1), 1, row_x, nbr_bottom, 7, &
                    f_post(7, 1, ny+1), 1, row_x, nbr_top, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, 1, 1), 1, row_x, nbr_bottom, 8, &
                    f_post(8, 1, ny+1), 1, row_x, nbr_top, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f_post(1, nx, 1), 1, column_y, nbr_right, 1, &
                    f_post(1, 0, 1), 1, column_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, nx, 1), 1, column_y, nbr_right, 5, &
                    f_post(5, 0, 1), 1, column_y, nbr_left, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, nx, 1), 1, column_y, nbr_right, 8, &
                    f_post(8, 0, 1), 1, column_y, nbr_left, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    ! message passing to left(i--)
    call MPI_Sendrecv(f_post(3, 1, 1), 1, column_y, nbr_left, 3, &
                    f_post(3, nx+1, 1), 1, column_y, nbr_right, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)       

    call MPI_Sendrecv(f_post(6, 1, 1), 1, column_y, nbr_left, 6, &
                    f_post(6, nx+1, 1), 1, column_y, nbr_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    call MPI_Sendrecv(f_post(7, 1, 1), 1, column_y, nbr_left, 7, &
                    f_post(7, nx+1, 1), 1, column_y, nbr_right, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Sendrecv(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, &
                    f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Sendrecv(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, &
                    f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Sendrecv(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, &
                    f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Sendrecv(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, &
                    f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

end subroutine message_passing_sendrecv