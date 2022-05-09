subroutine message_passing_nonblock()
    use mpi
    use commondata
    implicit none
    integer :: req(30)

    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(f_post(2, 1, ny), 1, row_x, nbr_top, 2, comm2d, req(1), rc)
    call MPI_Irecv(f_post(2, 1, 0), 1, row_x, nbr_bottom, 2, comm2d, req(2), rc)

    call MPI_Isend(f_post(5, 1, ny), 1, row_x, nbr_top, 5, comm2d, req(3), rc)
    call MPI_Irecv(f_post(5, 1, 0), 1, row_x, nbr_bottom, 5, comm2d, req(4), rc)

    call MPI_Isend(f_post(6, 1, ny), 1, row_x, nbr_top, 6, comm2d, req(5), rc)
    call MPI_Irecv(f_post(6, 1, 0), 1, row_x, nbr_bottom, 6, comm2d, req(6), rc)

    ! message passing to bottom(j--)
    call MPI_Isend(f_post(4, 1, 1), 1, row_x, nbr_bottom, 4, comm2d, req(7), rc)
    call MPI_Irecv(f_post(4, 1, ny+1), 1, row_x, nbr_top, 4, comm2d, req(8), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, row_x, nbr_bottom, 7, comm2d, req(9), rc)
    call MPI_Irecv(f_post(7, 1, ny+1), 1, row_x, nbr_top, 7, comm2d, req(10), rc)

    call MPI_Isend(f_post(8, 1, 1), 1, row_x, nbr_bottom, 8, comm2d, req(11), rc)
    call MPI_Irecv(f_post(8, 1, ny+1), 1, row_x, nbr_top, 8, comm2d, req(12), rc)


    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(f_post(1, nx, 1), 1, column_y, nbr_right, 1, comm2d, req(13), rc)
    call MPI_Irecv(f_post(1, 0, 1), 1, column_y, nbr_left, 1, comm2d, req(14), rc)

    call MPI_Isend(f_post(5, nx, 1), 1, column_y, nbr_right, 5, comm2d, req(15), rc)
    call MPI_Irecv(f_post(5, 0, 1), 1, column_y, nbr_left, 5, comm2d, req(16), rc)

    call MPI_Isend(f_post(8, nx, 1), 1, column_y, nbr_right, 8, comm2d, req(17), rc)
    call MPI_Irecv(f_post(8, 0, 1), 1, column_y, nbr_left, 8, comm2d, req(18), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(f_post(3, 1, 1), 1, column_y, nbr_left, 3, comm2d, req(19), rc)
    call MPI_Irecv(f_post(3, nx+1, 1), 1, column_y, nbr_right, 3, comm2d, req(20), rc)

    call MPI_Isend(f_post(6, 1, 1), 1, column_y, nbr_left, 6, comm2d, req(21), rc)
    call MPI_Irecv(f_post(6, nx+1, 1), 1, column_y, nbr_right, 6, comm2d, req(22), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, column_y, nbr_left, 7, comm2d, req(23), rc)
    call MPI_Irecv(f_post(7, nx+1, 1), 1, column_y, nbr_right, 7, comm2d, req(24), rc)   


    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Isend(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, comm2d, req(25), rc)
    call MPI_Irecv(f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, comm2d, req(26), rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Isend(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, comm2d, req(27), rc)
    call MPI_Irecv(f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, comm2d, req(28), rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Isend(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, comm2d, req(29), rc)
    call MPI_Irecv(f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, comm2d, req(30), rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Isend(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, comm2d, req(31), rc)
    call MPI_Irecv(f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, comm2d, req(32), rc)

    call MPI_Waitall(32, req, MPI_STATUSES_IGNORE, rc)
    
end subroutine message_passing_nonblock