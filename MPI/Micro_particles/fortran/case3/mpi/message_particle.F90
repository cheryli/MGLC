subroutine message_particle_f()
    use mpi
    use commondata
    implicit none
    integer :: cNum
    ! four corners, send to / recv from right(1), top(2), left(3), bottom(4)
    logical :: send_flag(4), recv_flag(4)  
    integer :: i, j, idx
    integer :: lf_b, rg_b, tp_b, bt_b
    integer :: start, len_x, len_y
    integer :: req(cNumMax * 64), cnt

    cnt = 0
    send_flag = .false.
    recv_flag = .false.

    ! assume effect zone is a square
    do cNum=1,cNumMax
        ! boundary of the zone
        lf_b = ceiling(xCenter(cNum) - radius(cNum) - 3) - i_start_global
        rg_b = floor(xCenter(cNum) + radius(cNum) + 3) - i_start_global
        bt_b = ceiling(yCenter(cNum) - radius(cNum) - 3) - j_start_global
        tp_b = floor(yCenter(cNum) + radius(cNum) + 3) - j_start_global

        ! particle effect this sub_domain
        if (-2 <= rg_b .AND. nx+3 >= lf_b .AND. &
            -2 <= tp_b .AND. ny+3 >= bt_b) then    

                start = max(bt_b, 1)
                len_y = max(min(tp_b, ny) - start + 1, 0)
                ! /*-------------------------- exchange message along x -----------------------------*/
                ! /*----------message passing to left(i--)
                ! /*------------- 1,2,3 -(to left)-> nx+1, nx+2, nx+3 ---------------*/
                ! left 1, 2, 3 line inside effect zone
                do i = 1, 3
                    if (i >= lf_b .AND. i <= rg_b) then
                        send_flag(3) = .true.
                        ! send 1,2,3 line to left neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, i, start), len_y, type_f_y, nbr_left, i, comm2d, req(cnt), rc)
                    endif
                enddo

                ! right nx+1, nx+2, nx+3 line inside effect zone, same tag as senders
                do i = 1, 3
                    if (nx+i >= lf_b .AND. nx+i <= rg_b) then
                        recv_flag(1) = .true.
                        ! receive right nx+1, nx+2, nx+3 lines from right neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, nx+i, start), len_y, type_f_y, nbr_right, i, comm2d, req(cnt), rc)
                    endif
                enddo

                ! /*----------message passing to right(i++)
                ! /*------------- nx, nx-1, nx-2 -(to right)-> 0, -1, -2 ---------------*/
                ! right nx, nx-1, nx-2 line inside effect zone
                do i = 1, 3
                    if (nx+1-i >= lf_b .AND. nx+1-i <= rg_b) then
                        send_flag(1) = .true.
                        ! send nx, nx-1, nx-2 line to right neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, nx+1-i, start), len_y, type_f_y, nbr_right, 10+i, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! left 0, -1, -2 line inside effect zone, same tag as senders
                do i = 1, 3
                    if (1-i >= lf_b .AND. 1-i <= rg_b) then
                        recv_flag(3) = .true.
                        ! receive left 0, -1, -2 lines from left neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, 1-i, start), len_y, type_f_y, nbr_left, 10+i, comm2d, req(cnt), rc)
                    endif
                enddo

                start = max(lf_b, 1)
                len_x = max(min(rg_b, nx) - start + 1, 0)
                ! /*-------------------------- exchange message along y -----------------------------*/
                ! /*----------message passing to bottom(j--)
                ! /*------------------- 1,2,3 -(to bottom)-> ny+1, ny+2, ny+3 ---------------------*/
                ! bottom line 1, 2, 3 inside effect zone
                do j = 1, 3
                    if (j >= bt_b .AND. j <= tp_b) then
                        send_flag(4) = .true.
                        ! send 1,2,3 line to bottom neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, start, j), len_x, type_f_x, nbr_bottom, 20+j, comm2d, req(cnt), rc)
                    endif
                enddo

                ! top ny+1, ny+2, ny+3 line inside effect zone, same tag as senders
                do j = 1, 3
                    if (ny+j >= bt_b .AND. ny+j <= tp_b) then
                        recv_flag(2) = .true.
                        ! receive top ny+1, ny+2, ny+3 lines from top neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, start, ny+j), len_x, type_f_x, nbr_top, 20+j, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! /*----------message passing to top(j++)
                ! /*------------------- ny, ny-1, ny-2 -(to top)-> 0, -1, -2 ---------------------*/
                ! top ny, ny-1, ny-2 line inside effect zone
                do j = 1, 3
                    if (ny+1-j >= bt_b .AND. ny+1-j <= tp_b) then
                        send_flag(2) = .true.
                        ! send ny, ny-1, ny-2 line to top neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, start, ny+1-j), len_x, type_f_x, nbr_top, 30+j, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! bottom 0, -1, -2 line inside effect zone, same tag as senders
                do j = 1, 3
                    if (1-j >= bt_b .AND. 1-j <= tp_b) then
                        recv_flag(4) = .true.
                        ! receive bottom 0, -1, -2 lines from bottom neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, start, 1-j), len_x, type_f_x, nbr_bottom, 30+j, comm2d, req(cnt), rc)
                    endif
                enddo



                ! /*-------------------------- corner messages -----------------------------*/
                ! using pipline to send multiple times or send corners as a whole
                ! message passing to top-right corner(i++, j++)
                ! send to top-right corner (1, 2)
                if (send_flag(1) .AND. send_flag(2)) then
                    cnt = cnt + 1
                    call MPI_Isend(f(0, nx-2, ny-2), 1, type_square3_f, cnr_top_right, 41, comm2d, req(cnt), rc)
                endif
                ! receive from bottom-left corner (3, 4)
                if (recv_flag(3) .AND. recv_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f(0, -2, -2), 1, type_square3_f, cnr_bottom_left, 41, comm2d, req(cnt), rc)
                endif

                ! message passing to top-left corner(i--, j++)
                ! send to top-left corner (2, 3)
                if (send_flag(2) .AND. send_flag(3)) then
                    cnt = cnt + 1
                    call MPI_Isend(f(0, 1, ny-2), 1, type_square3_f, cnr_top_left, 42, comm2d, req(cnt), rc)
                endif
                ! receive from bottom-right corner (1, 4)
                if (recv_flag(1) .AND. recv_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f(0, nx+1, -2), 1, type_square3_f, cnr_bottom_right, 42, comm2d, req(cnt), rc)
                endif
                
                ! message passing to bottom-left corner(i--, j--)
                ! send to bottom-left corner (3, 4)
                if (send_flag(3) .AND. send_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Isend(f(0, 1, 1), 1, type_square3_f, cnr_bottom_left, 43, comm2d, req(cnt), rc)
                endif
                ! receive from top-right corner (1, 2)
                if (recv_flag(1) .AND. recv_flag(2)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f(0, nx+1, ny+1), 1, type_square3_f, cnr_top_right, 43, comm2d, req(cnt), rc)
                endif

                ! message passing to bottom-right corner(i++, j--)
                ! send to bottom-right corner (1, 4)
                if (send_flag(1) .AND. send_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Isend(f(0, nx-2, 1), 1, type_square3_f, cnr_bottom_right, 44, comm2d, req(cnt), rc)
                endif
                ! receive from top-left corner (2, 3)
                if (recv_flag(2) .AND. recv_flag(3)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f(0, -2, ny+1), 1, type_square3_f, cnr_top_left, 44, comm2d, req(cnt), rc)
                endif

        endif
    enddo

    call MPI_Waitall(cnt, req, MPI_STATUSES_IGNORE, rc)

end subroutine message_particle_f




subroutine message_particle_f_post()
    use mpi
    use commondata
    implicit none
    integer :: cNum
    ! four corners, send to / recv from right(1), top(2), left(3), bottom(4)
    logical :: send_flag(4), recv_flag(4)
    integer :: i
    integer :: lf_b, rg_b, tp_b, bt_b
    integer :: start, len_x, len_y
    integer :: req(cNumMax * 64), cnt

    cnt = 0
    send_flag = .false.
    recv_flag = .false.

    ! assume effect zone is a square
    do cNum=1,cNumMax
        ! boundary of the zone
        lf_b = ceiling(xCenter(cNum) - radius(cNum) - 3) - i_start_global
        rg_b = floor(xCenter(cNum) + radius(cNum) + 3) - i_start_global
        bt_b = ceiling(yCenter(cNum) - radius(cNum) - 3) - j_start_global
        tp_b = floor(yCenter(cNum) + radius(cNum) + 3) - j_start_global

        ! particle effect this sub_domain
        if (-1 <= rg_b .AND. nx+2 >= lf_b .AND. &
            -1 <= tp_b .AND. ny+2 >= bt_b) then    

                start = max(bt_b, 1)
                len_y = max(min(tp_b, ny) - start + 1, 0)
                ! /*-------------------------- exchange message along x -----------------------------*/
                ! /*----------message passing to left(i--)
                ! /*------------------- 1 -(to left)-> nx+1 ---------------------*/
                ! left 1 line inside effect zone
                if (1 >= lf_b .AND. 1 <= rg_b) then
                    send_flag(3) = .true.
                    ! send fp_1,5,8 (1) to left neighbor
                    call MPI_Isend(f_post(1, 1, start), len_y, type_fpi_y, nbr_left, 1, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, 1, start), len_y, type_fpi_y, nbr_left, 5, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, 1, start), len_y, type_fpi_y, nbr_left, 8, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! right nx+1(neibor 1) line inside effect zone, same tag as senders
                if (nx+1 >= lf_b .AND. nx+1 <= rg_b) then
                    recv_flag(1) = .true.
                    ! receive fp_1,5,8 (nx+1) from right neighbor 
                    call MPI_Irecv(f_post(1, nx+1, start), len_y, type_fpi_y, nbr_right, 1, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, nx+1, start), len_y, type_fpi_y, nbr_right, 5, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, nx+1, start), len_y, type_fpi_y, nbr_right, 8, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- 2 -(to left)-> nx+2 ---------------------*/
                ! left 2 line inside effect zone
                if (2 >= lf_b .AND. 2 <= rg_b) then
                    send_flag(3) = .true.
                    ! send fp_3,6,7 (2) to left neighbor
                    call MPI_Isend(f_post(3, 2, start), len_y, type_fpi_y, nbr_left, 3, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(6, 2, start), len_y, type_fpi_y, nbr_left, 6, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(7, 2, start), len_y, type_fpi_y, nbr_left, 7, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! right nx+2(neibor 2) line inside effect zone, same tag as senders
                if (nx+2 >= lf_b .AND. nx+2 <= rg_b) then
                    recv_flag(1) = .true.
                    ! receive fp_3,6,7 (nx+2) from right neighbor 
                    call MPI_Irecv(f_post(3, nx+2, start), len_y, type_fpi_y, nbr_right, 3, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(6, nx+2, start), len_y, type_fpi_y, nbr_right, 6, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(7, nx+2, start), len_y, type_fpi_y, nbr_right, 7, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*----------message passing to right(i++)
                ! /*------------------- nx -(to right)-> 0 ---------------------*/
                ! right nx line inside effect zone
                if (nx >= lf_b .AND. nx <= rg_b) then
                    send_flag(1) = .true.
                    ! send fp_3,6,7 (nx) to right neighbor
                    call MPI_Isend(f_post(3, nx, start), len_y, type_fpi_y, nbr_right, 13, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(6, nx, start), len_y, type_fpi_y, nbr_right, 16, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(7, nx, start), len_y, type_fpi_y, nbr_right, 17, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! left 0 line inside effect zone
                if (0 >= lf_b .AND. 0 <= rg_b) then
                    recv_flag(3) = .true.
                    ! receive fp_3,6,7 (0) from left neighbor
                    call MPI_Irecv(f_post(3, 0, start), len_y, type_fpi_y, nbr_left, 13, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(6, 0, start), len_y, type_fpi_y, nbr_left, 16, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(7, 0, start), len_y, type_fpi_y, nbr_left, 17, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- nx-1 -(to right)-> -1 ---------------------*/
                if (nx-1 >= lf_b .AND. nx-1 <= rg_b) then
                    send_flag(1) = .true.
                    ! send fp_1,5,8 (nx-1) to right neighbor
                    call MPI_Isend(f_post(1, nx-1, start), len_y, type_fpi_y, nbr_right, 11, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, nx-1, start), len_y, type_fpi_y, nbr_right, 15, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, nx-1, start), len_y, type_fpi_y, nbr_right, 18, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! left -1 line inside effect zone
                if (-1 >= lf_b .AND. -1 <= rg_b) then
                    recv_flag(3) = .true.
                    ! receive fp_1,5,8 (-1) from left neighbor
                    call MPI_Irecv(f_post(1, -1, start), len_y, type_fpi_y, nbr_left, 11, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, -1, start), len_y, type_fpi_y, nbr_left, 15, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, -1, start), len_y, type_fpi_y, nbr_left, 18, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif


                start = max(lf_b, 1)
                len_x = max(min(rg_b, nx) - start + 1, 0)
                ! /*-------------------------- exchange message along y -----------------------------*/
                ! /*----------message passing to bottom(j--)
                ! /*------------------- 1 -(to bottom)-> ny+1 ---------------------*/
                ! bottom 1 line inside effect zone
                if (1 >= bt_b .AND. 1 <= tp_b) then
                    send_flag(4) = .true.
                    ! send fp_2,5,6 (1) to bottom neighbor
                    call MPI_Isend(f_post(2, start, 1), len_x, type_fi_x, nbr_bottom, 22, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, start, 1), len_x, type_fi_x, nbr_bottom, 25, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(6, start, 1), len_x, type_fi_x, nbr_bottom, 26, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! top ny+1(neibor 1) line inside effect zone, same tag as senders
                if (ny+1 >= bt_b .AND. ny+1 <= tp_b) then
                    recv_flag(2) = .true.
                    ! receive fp_2,5,6 (ny+1) lines from top neighbor 
                    call MPI_Irecv(f_post(2, start, ny+1), len_x, type_fi_x, nbr_top, 22, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, start, ny+1), len_x, type_fi_x, nbr_top, 25, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(6, start, ny+1), len_x, type_fi_x, nbr_top, 26, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- 2 -(to bottom)-> ny+2 ---------------------*/
                ! bottom 2 line inside effect zone
                if (2 >= bt_b .AND. 2 <= tp_b) then
                    send_flag(4) = .true.
                    ! send fp_4,7,8 (2) to bottom neighbor
                    call MPI_Isend(f_post(4, start, 2), len_x, type_fi_x, nbr_bottom, 24, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(7, start, 2), len_x, type_fi_x, nbr_bottom, 27, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, start, 2), len_x, type_fi_x, nbr_bottom, 28, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! top ny+2(neibor 2) line inside effect zone, same tag as senders
                if (ny+2 >= bt_b .AND. ny+2 <= tp_b) then
                    recv_flag(2) = .true.
                    ! receive fp_4,7,8 (ny+1) lines from top neighbor 
                    call MPI_Irecv(f_post(4, start, ny+2), len_x, type_fi_x, nbr_top, 24, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(7, start, ny+2), len_x, type_fi_x, nbr_top, 27, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, start, ny+2), len_x, type_fi_x, nbr_top, 28, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif


                ! /*----------message passing to top(j++)
                ! /*------------------- ny -(to top)-> 0 ---------------------*/
                ! top ny line inside effect zone
                if (ny >= bt_b .AND. ny <= tp_b) then
                    send_flag(2) = .true.
                    ! send fp_4,7,8 (ny) to top neighbor
                    call MPI_Isend(f_post(4, start, ny), len_x, type_fi_x, nbr_top, 34, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(7, start, ny), len_x, type_fi_x, nbr_top, 37, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, start, ny), len_x, type_fi_x, nbr_top, 38, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! bottom 0 line inside effect zone
                if (0 >= bt_b .AND. 0 <= tp_b) then
                    recv_flag(4) = .true.
                    ! receive fp_4,7,8 (0) from bottom neighbor
                    call MPI_Irecv(f_post(4, start, 0), len_x, type_fi_x, nbr_bottom, 34, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(7, start, 0), len_x, type_fi_x, nbr_bottom, 37, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, start, 0), len_x, type_fi_x, nbr_bottom, 38, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                ! /*------------------- ny-1 -(to top)-> -1 ---------------------*/
                ! top ny-1 line inside effect zone
                if (ny-1 >= bt_b .AND. ny-1 <= tp_b) then
                    send_flag(2) = .true.
                    ! send fp_2,5,6 (ny-1) to top neighbor
                    call MPI_Isend(f_post(2, start, ny-1), len_x, type_fi_x, nbr_top, 32, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, start, ny-1), len_x, type_fi_x, nbr_top, 35, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(6, start, ny-1), len_x, type_fi_x, nbr_top, 36, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! bottom -1 line inside effect zone
                if (-1 >= bt_b .AND. -1 <= tp_b) then
                    recv_flag(4) = .true.
                    ! receive fp_2,5,6 (-1) from bottom neighbor
                    call MPI_Irecv(f_post(2, start, -1), len_x, type_fi_x, nbr_bottom, 32, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, start, -1), len_x, type_fi_x, nbr_bottom, 35, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(6, start, -1), len_x, type_fi_x, nbr_bottom, 36, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif


                ! /*-------------------------- corner messages -----------------------------*/
                ! using pipline to send multiple times or send corners as a whole
                ! message passing to top-right corner(i++, j++)
                ! send to top-right corner (1, 2)
                if (send_flag(1) .AND. send_flag(2)) then
                    cnt = cnt + 1
                    call MPI_Isend(f_post(7, nx, ny), 1, MPI_REAL8, cnr_top_right, 41, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Isend(f_post(5, nx-1, ny-1), 1, type_square2_fpi, cnr_top_right, 42, comm2d, req(cnt), rc)
                endif
                ! receive from bottom-left corner (3, 4)
                if (recv_flag(3) .AND. recv_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(7, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 41, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(5, -1, -1), 1, type_square2_fpi, cnr_bottom_left, 42, comm2d, req(cnt), rc)
                endif

                ! message passing to top-left corner(i--, j++)
                ! send to top-left corner (2, 3)
                if (send_flag(2) .AND. send_flag(3)) then
                    cnt = cnt + 1
                    call MPI_Isend(f_post(8, 1, ny), 1, MPI_REAL8, cnr_top_left, 43, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Isend(f_post(6, 1, ny-1), 1, type_square2_fpi, cnr_top_left, 44, comm2d, req(cnt), rc)
                endif
                ! receive from bottom-right corner (1, 4)
                if (recv_flag(1) .AND. recv_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(8, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 43, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(6, nx+1, -1), 1, type_square2_fpi, cnr_bottom_right, 44, comm2d, req(cnt), rc)
                endif
            
                ! message passing to bottom-left corner(i--, j--)
                ! send to bottom-left corner (3, 4)
                if (send_flag(3) .AND. send_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Isend(f_post(5, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 45, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Isend(f_post(7, 1, 1), 1, type_square2_fpi, cnr_bottom_left, 46, comm2d, req(cnt), rc)
                endif
                ! receive from top-right corner (1, 2)
                if (recv_flag(1) .AND. recv_flag(2)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(5, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 45, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(7, nx+1, ny+1), 1, type_square2_fpi, cnr_top_right, 46, comm2d, req(cnt), rc)
                endif

                ! message passing to bottom-right corner(i++, j--)
                ! send to bottom-right corner (1, 4)
                if (send_flag(1) .AND. send_flag(4)) then
                    cnt = cnt + 1
                    call MPI_Isend(f_post(6, nx, 1), 1, type_square2_fpi, cnr_bottom_right, 47, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Isend(f_post(8, nx-1, 1), 1, type_square2_fpi, cnr_bottom_right, 48, comm2d, req(cnt), rc)
                endif
                ! receive from top-left corner (2, 3)
                if (recv_flag(2) .AND. recv_flag(3)) then
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(6, 0, ny+1), 1, type_square2_fpi, cnr_top_left, 47, comm2d, req(cnt), rc)
                    cnt = cnt + 1
                    call MPI_Irecv(f_post(8, -1, ny+1), 1, type_square2_fpi, cnr_top_left, 48, comm2d, req(cnt), rc)
                endif

                ! ! /*-------------------------- corner messages -----------------------------*/
                ! sending all fp
                ! ! using pipline to send multiple times or send corners as a whole
                ! ! message passing to top-right corner(i++, j++)
                ! ! send to top-right corner (1, 2)
                ! if (send_flag(1) .AND. send_flag(2)) then
                !     cnt = cnt + 1
                !     call MPI_Isend(f_post(0, nx-1, ny-1), 1, type_square2_fp, cnr_top_right, 41, comm2d, req(cnt), rc)
                ! endif
                ! ! receive from bottom-left corner (3, 4)
                ! if (recv_flag(3) .AND. recv_flag(4)) then
                !     cnt = cnt + 1
                !     call MPI_Irecv(f_post(0, -1, -1), 1, type_square2_fp, cnr_bottom_left, 41, comm2d, req(cnt), rc)
                ! endif

                ! ! message passing to top-left corner(i--, j++)
                ! ! send to top-left corner (2, 3)
                ! if (send_flag(2) .AND. send_flag(3)) then
                !     cnt = cnt + 1
                !     call MPI_Isend(f_post(0, 1, ny-1), 1, type_square2_fp, cnr_top_left, 42, comm2d, req(cnt), rc)
                ! endif
                ! ! receive from bottom-right corner (1, 4)
                ! if (recv_flag(1) .AND. recv_flag(4)) then
                !     cnt = cnt + 1
                !     call MPI_Irecv(f_post(0, nx+1, -1), 1, type_square2_fp, cnr_bottom_right, 42, comm2d, req(cnt), rc)
                ! endif
            
                ! ! message passing to bottom-left corner(i--, j--)
                ! ! send to bottom-left corner (3, 4)
                ! if (send_flag(3) .AND. send_flag(4)) then
                !     cnt = cnt + 1
                !     call MPI_Isend(f_post(0, 1, 1), 1, type_square2_fp, cnr_bottom_left, 43, comm2d, req(cnt), rc)
                ! endif
                ! ! receive from top-right corner (1, 2)
                ! if (recv_flag(1) .AND. recv_flag(2)) then
                !     cnt = cnt + 1
                !     call MPI_Irecv(f_post(0, nx+1, ny+1), 1, type_square2_fp, cnr_top_right, 43, comm2d, req(cnt), rc)
                ! endif

                ! ! message passing to bottom-right corner(i++, j--)
                ! ! send to bottom-right corner (1, 4)
                ! if (send_flag(1) .AND. send_flag(4)) then
                !     cnt = cnt + 1
                !     call MPI_Isend(f_post(0, nx-1, 1), 1, type_square2_fp, cnr_bottom_right, 44, comm2d, req(cnt), rc)
                ! endif
                ! ! receive from top-left corner (2, 3)
                ! if (recv_flag(2) .AND. recv_flag(3)) then
                !     cnt = cnt + 1
                !     call MPI_Irecv(f_post(0, -1, ny+1), 1, type_square2_fp, cnr_top_left, 44, comm2d, req(cnt), rc)
                ! endif
                
        endif
    enddo

    call MPI_Waitall(cnt, req, MPI_STATUSES_IGNORE, rc)

end subroutine message_particle_f_post