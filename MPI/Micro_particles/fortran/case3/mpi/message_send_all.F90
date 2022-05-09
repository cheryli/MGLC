subroutine send_all_f()
    use mpi
    use commondata
    implicit none

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f(0, nx, 1), ny, type_f_y, nbr_right, 1, &
                    f(0, 0, 1), ny, type_f_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, nx-1, 1), ny, type_f_y, nbr_right, 2, &
                    f(0, -1, 1), ny, type_f_y, nbr_left, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, nx-2, 1), ny, type_f_y, nbr_right, 3, &
                    f(0, -2, 1), ny, type_f_y, nbr_left, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to left(i--) 
    call MPI_Sendrecv(f(0, 1, 1), ny, type_f_y, nbr_left, 4, &
                    f(0, nx+1, 1), ny, type_f_y, nbr_right, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)          
                
    call MPI_Sendrecv(f(0, 2, 1), ny, type_f_y, nbr_left, 5, &
                    f(0, nx+2, 1), ny, type_f_y, nbr_right, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)   

    call MPI_Sendrecv(f(0, 3, 1), ny, type_f_y, nbr_left, 6, &
                    f(0, nx+3, 1), ny, type_f_y, nbr_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)   

    ! message passing to top(j++)
    call MPI_Sendrecv(f(0, 1, ny), nx, type_f_x, nbr_top, 7, &
                    f(0, 1, 0), nx, type_f_x, nbr_bottom, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    
    call MPI_Sendrecv(f(0, 1, ny-1), nx, type_f_x, nbr_top, 8, &
                    f(0, 1, -1), nx, type_f_x, nbr_bottom, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, 1, ny-2), nx, type_f_x, nbr_top, 9, &
                    f(0, 1, -2), nx, type_f_x, nbr_bottom, 9, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to bottom(j--)
    call MPI_Sendrecv(f(0, 1, 1), nx, type_f_x, nbr_bottom, 10, &
                    f(0, 1, ny+1), nx, type_f_x, nbr_top, 10, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, 1, 2), nx, type_f_x, nbr_bottom, 11, &
                    f(0, 1, ny+2), nx, type_f_x, nbr_top, 11, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, 1, 3), nx, type_f_x, nbr_bottom, 12, &
                    f(0, 1, ny+3), nx, type_f_x, nbr_top, 12, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! corner messages
    ! exchange message with corners --- message tag: discrete velocity
    ! ! exchange message along top-right (i++, j++)
    call MPI_Sendrecv(f(0, nx-2, ny-2), 1, type_square3_f, cnr_top_right, 15, &
                    f(0, -2, -2), 1, type_square3_f, cnr_bottom_left, 15, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ! exchange message along top-left (i--, j++)
    call MPI_Sendrecv(f(0, 1, ny-2), 1, type_square3_f, cnr_top_left, 16, &
                    f(0, nx+1, -2), 1, type_square3_f, cnr_bottom_right, 16, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ! exchange message along bottom-left (i--, j--)
    call MPI_Sendrecv(f(0, 1, 1), 1, type_square3_f, cnr_bottom_left, 17, &
                    f(0, nx+1, ny+1), 1, type_square3_f, cnr_top_right, 17, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ! exchange message along bottom-right (i++, j--)
    call MPI_Sendrecv(f(0, nx-2, 1), 1, type_square3_f, cnr_bottom_right, 18, &
                f(0, -2, ny+1), 1, type_square3_f, cnr_top_left, 18, &
                comm2d, MPI_STATUS_IGNORE, rc)

end subroutine send_all_f


subroutine send_all_fp()
    use mpi
    use commondata
    implicit none

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f_post(0, nx, 1), ny, type_fp_y, nbr_right, 1, &
                    f_post(0, 0, 1), ny, type_fp_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(0, nx-1, 1), ny, type_fp_y, nbr_right, 2, &
                    f_post(0, -1, 1), ny, type_fp_y, nbr_left, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to left(i--)
    call MPI_Sendrecv(f_post(0, 1, 1), ny, type_fp_y, nbr_left, 3, &
                    f_post(0, nx+1, 1), ny, type_fp_y, nbr_right, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)    

    call MPI_Sendrecv(f_post(0, 2, 1), ny, type_fp_y, nbr_left, 4, &
                    f_post(0, nx+2, 1), ny, type_fp_y, nbr_right, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    ! message passing to top(j++)
    call MPI_Sendrecv(f_post(0, 1, ny), nx, type_f_x, nbr_top, 5, &
                    f_post(0, 1, 0), nx, type_f_x, nbr_bottom, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    
    call MPI_Sendrecv(f_post(0, 1, ny-1), nx, type_f_x, nbr_top, 6, &
                    f_post(0, 1, -1), nx, type_f_x, nbr_bottom, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    
    ! message passing to bottom(j--)
    call MPI_Sendrecv(f_post(0, 1, 1), nx, type_f_x, nbr_bottom, 7, &
                    f_post(0, 1, ny+1), nx, type_f_x, nbr_top, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    
    call MPI_Sendrecv(f_post(0, 1, 2), nx, type_f_x, nbr_bottom, 8, &
                    f_post(0, 1, ny+2), nx, type_f_x, nbr_top, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)


    ! exchange message with corners
    ! exchange message along top-right (i++, j++)
    call MPI_Sendrecv(f_post(0, nx-1, ny-1), 1, type_square2_fp, cnr_top_right, 15, &
                    f_post(0, -1, -1), 1, type_square2_fp, cnr_bottom_left, 15, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Sendrecv(f_post(0, 1, ny-1), 1, type_square2_fp, cnr_top_left, 16, &
                    f_post(0, nx+1, -1), 1, type_square2_fp, cnr_bottom_right, 16, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Sendrecv(f_post(0, 1, 1), 1, type_square2_fp, cnr_bottom_left, 17, &
                    f_post(0, nx+1, ny+1), 1, type_square2_fp, cnr_top_right, 17, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Sendrecv(f_post(0, nx-1, 1), 1, type_square2_fp, cnr_bottom_right, 18, &
                f_post(0, -1, ny+1), 1, type_square2_fp, cnr_top_left, 18, &
                comm2d, MPI_STATUS_IGNORE, rc)

end subroutine send_all_fp