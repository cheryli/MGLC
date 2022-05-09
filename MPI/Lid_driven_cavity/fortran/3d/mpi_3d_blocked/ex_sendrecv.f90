subroutine message_passing_sendrecv()
    use mpi
    use commondata
    implicit none
    integer :: tag(5), idx

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)s
    tag = (/ 1, 7, 9, 11, 13 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), nx, 1, 1), 1, surface_x, nbr_surface(1), tag(idx), &
                    f_post(tag(idx), 0, 1, 1), 1, surface_x, nbr_surface(2), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (i--)
    tag = (/ 2, 8, 10 ,12, 14 /) 
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_x, nbr_surface(2), tag(idx), &
                    f_post(tag(idx), nx+1, 1, 1), 1, surface_x, nbr_surface(1), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    tag = (/ 3, 7, 8, 15, 17 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, ny, 1), 1, surface_y, nbr_surface(3), tag(idx), &
                    f_post(tag(idx), 1, 0, 1), 1, surface_y, nbr_surface(4), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (j--)
    tag = (/ 4, 9, 10, 16, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_y, nbr_surface(4), tag(idx), &
                    f_post(tag(idx), 1, ny+1, 1), 1, surface_y, nbr_surface(3), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    tag = (/ 5, 11, 12, 15, 16 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, nz), 1, surface_z, nbr_surface(5), tag(idx), &
                    f_post(tag(idx), 1, 1, 0), 1, surface_z, nbr_surface(6), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (k--)
    tag = (/ 6, 13, 14, 17, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_z, nbr_surface(6), tag(idx), &
                    f_post(tag(idx), 1, 1, nz+1), 1, surface_z, nbr_surface(5), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! /* ----------------------------- line data ----------------------------- */
    ! ------------ exchange message perpendicular to z ----------------
    ! message passing to (i++, j++) (7)
    call MPI_Sendrecv(f_post(7, nx, ny, 1), 1, line_z, nbr_line(7), 7, &
                    f_post(7, 0, 0, 1), 1, line_z, nbr_line(10), 7, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j--) (10)
    call MPI_Sendrecv(f_post(10, 1, 1, 1), 1, line_z, nbr_line(10), 10, &
                    f_post(10, nx+1, ny+1, 1), 1, line_z, nbr_line(7), 10, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, j--) (9)
    call MPI_Sendrecv(f_post(9, nx, 1, 1), 1, line_z, nbr_line(9), 9, &
                    f_post(9, 0, ny+1, 1), 1, line_z, nbr_line(8), 9, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j++) (8)
    call MPI_Sendrecv(f_post(8, 1, ny, 1), 1, line_z, nbr_line(8), 8, &
                    f_post(8, nx+1, 0, 1), 1, line_z, nbr_line(9), 8, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to y ----------------
    ! message passing to (i++, k++) (11)
    call MPI_Sendrecv(f_post(11, nx, 1, nz), 1, line_y, nbr_line(11), 11, &
                    f_post(11, 0, 1, 0), 1, line_y, nbr_line(14), 11, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k--) (14)
    call MPI_Sendrecv(f_post(14, 1, 1, 1), 1, line_y, nbr_line(14), 14, &
                    f_post(14, nx+1, 1, nz+1), 1, line_y, nbr_line(11), 14, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, k--) (13)
    call MPI_Sendrecv(f_post(13, nx, 1, 1), 1, line_y, nbr_line(13), 13, &
                    f_post(13, 0, 1, nz+1), 1, line_y, nbr_line(12), 13, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k++) (12)
    call MPI_Sendrecv(f_post(12, 1, 1, nz), 1, line_y, nbr_line(12), 12, &
                    f_post(12, nx+1, 1, 0), 1, line_y, nbr_line(13), 12, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to x ----------------
    ! message passing to (j++, k++) (15)
    call MPI_Sendrecv(f_post(15, 1, ny, nz), 1, line_x, nbr_line(15), 15, &
                    f_post(15, 1, 0, 0), 1, line_x, nbr_line(18), 15, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k--) (18)
    call MPI_Sendrecv(f_post(18, 1, 1, 1), 1, line_x, nbr_line(18), 18, &
                    f_post(18, 1, ny+1, nz+1), 1, line_x, nbr_line(15), 18, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j++, k--) (17)
    call MPI_Sendrecv(f_post(17, 1, ny, 1), 1, line_x, nbr_line(17), 17, &
                    f_post(17, 1, 0, nz+1), 1, line_x, nbr_line(16), 17, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k++) (16)
    call MPI_Sendrecv(f_post(16, 1, 1, nz), 1, line_x, nbr_line(16), 16, &
                    f_post(16, 1, ny+1, 0), 1, line_x, nbr_line(17), 16, &
                    comm3d, MPI_STATUS_IGNORE, rc)
    

    ! don't need to echange point data


end subroutine message_passing_sendrecv