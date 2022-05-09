subroutine message_particle_f()
    use mpi
    use commondata
    implicit none
    integer :: cNum
    integer :: send_flag(4)  ! four corners
    integer :: i, j, idx_global
    integer :: lf_b, rg_b, tp_b, bt_b, len_x, len_y
    integer :: req(cNumMax * 64), cnt

    cnt = 0

    ! assume effect zone is a square
    do cNum=1,cNumMax
        ! boundary of the zone
        lf_b = ceiling(xCenter(cNum) - radius(cNum) - 3)
        rg_b = floor(xCenter(cNum) + radius(cNum) + 3)
        bt_b = ceiling(yCenter(cNum) - radius(cNum) - 3)
        tp_b = floor(yCenter(cNum) + radius(cNum) + 3)
        len_x = rg_b - lf_b + 1
        len_y = tp_b - bt_b + 1

        ! particle effect this sub_domain
        if (i_start_global - 2 <= rg_b .AND. &
            i_start_global + nx + 3 >= lf_b .AND. &
            j_start_global - 2 <= tp_b .AND. &
            j_start_global + ny + 3 >= bt_b) then    

                ! /*-------------------------- exchange message along x -----------------------------*/
                ! /*----------message passing to left(i--)
                ! /*------------- 1,2,3 -(to left)-> nx+1, nx+2, nx+3 ---------------*/
                ! left 1, 2, 3 line inside effect zone
                do i = 1, 3
                    idx_global = i_start_global + i
                    if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                        ! send 1,2,3 line to left neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, i, bt_b), len_y, type_f_y, nbr_left, i, comm2d, req(cnt), rc)
                    endif
                enddo

                ! right nx+1, nx+2, nx+3 line inside effect zone, same tag as senders
                do i = 1, 3
                    idx_global = i_start_global + nx + i
                    if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                        ! receive right nx+1, nx+2, nx+3 lines from right neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, nx+i, bt_b), len_y, type_f_y, nbr_right, i, comm2d, req(cnt), rc)
                    endif
                enddo

                ! /*----------message passing to right(i++)
                ! /*------------- nx, nx-1, nx-2 -(to right)-> 0, -1, -2 ---------------*/
                ! right nx, nx-1, nx-2 line inside effect zone
                do i = 1, 3
                    idx_global = i_start_global + nx + 1 - i
                    if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                        ! send nx, nx-1, nx-2 line to right neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, nx+1-i, bt_b), len_y, type_f_y, nbr_right, 10+i, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! left 0, -1, -2 line inside effect zone, same tag as senders
                do i = 1, 3
                    idx_global = i_start_global + 1 - i
                    if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                        ! receive left 0, -1, -2 lines from left neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, 1-i, bt_b), len_y, type_f_y, nbr_left, 10+i, comm2d, req(cnt), rc)
                    endif
                enddo


                ! /*-------------------------- exchange message along y -----------------------------*/
                ! /*----------message passing to bottom(j--)
                ! /*------------------- 1,2,3 -(to bottom)-> ny+1, ny+2, ny+3 ---------------------*/
                ! bottom line 1, 2, 3 inside effect zone
                do j = 1, 3
                    idx_global = j_start_global + j
                    if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                        ! send 1,2,3 line to bottom neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, lf_b, j), len_x, type_f_x, nbr_bottom, 20+j, comm2d, req(cnt), rc)
                    endif
                enddo

                ! top ny+1, ny+2, ny+3 line inside effect zone, same tag as senders
                do j = 1, 3
                    idx_global = j_start_global + ny + j
                    if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                        ! receive top ny+1, ny+2, ny+3 lines from top neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, lf_b, ny+j), len_x, type_f_x, nbr_top, 20+j, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! /*----------message passing to top(j++)
                ! /*------------------- ny, ny-1, ny-2 -(to top)-> 0, -1, -2 ---------------------*/
                ! top ny, ny-1, ny-2 line inside effect zone
                do j = 1, 3
                    idx_global = j_start_global + ny + 1 - j
                    if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                        ! send ny, ny-1, ny-2 line to right neighbor
                        cnt = cnt + 1
                        call MPI_Isend(f(0, lf_b, ny+1-j), len_x, type_f_x, nbr_top, 30+j, comm2d, req(cnt), rc)
                    endif
                enddo
                
                ! bottom 0, -1, -2 line inside effect zone, same tag as senders
                do j = 1, 3
                    idx_global = j_start_global + 1 - j
                    if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                        ! receive bottom 0, -1, -2 lines from left neighbor
                        cnt = cnt + 1
                        call MPI_Irecv(f(0, lf_b, 1-j), len_x, type_f_x, nbr_bottom, 30+j, comm2d, req(cnt), rc)
                    endif
                enddo



                ! /*-------------------------- corner messages -----------------------------*/
                
            
        endif
    enddo

    call MPI_Waitall(cnt, req, MPI_STATUSES_IGNORE, rc)

end subroutine message_particle_f




subroutine message_particle_f_post()
    use mpi
    use commondata
    implicit none
    integer :: cNum
    integer :: flag(4)  ! four corners
    integer :: i, idx_global
    integer :: lf_b, rg_b, tp_b, bt_b, len_x, len_y
    integer :: req(cNumMax * 64), cnt

    cnt = 0

    ! assume effect zone is a square
    do cNum=1,cNumMax
        ! boundary of the zone
        lf_b = ceiling(xCenter(cNum) - radius(cNum) - 3)
        rg_b = floor(xCenter(cNum) + radius(cNum) + 3)
        bt_b = ceiling(yCenter(cNum) - radius(cNum) - 3)
        tp_b = floor(yCenter(cNum) + radius(cNum) + 3)
        len_x = rg_b - lf_b + 1
        len_y = tp_b - bt_b + 1

        ! particle effect this sub_domain
        if (i_start_global - 1 <= rg_b .AND. &
            i_start_global + nx + 2 >= lf_b .AND. &
            j_start_global - 1 <= tp_b .AND. &
            j_start_global + ny + 2 >= bt_b) then    
                
                ! /*-------------------------- exchange message along x -----------------------------*/
                ! /*----------message passing to left(i--)
                ! /*------------------- 1 -(to left)-> nx+1 ---------------------*/
                ! left 1 line inside effect zone
                idx_global = i_start_global + 1
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! send fp_1,5,8 (1) to left neighbor
                    call MPI_Isend(f_post(1, 1, bt_b), len_y, type_fpi_y, nbr_left, 1, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, 1, bt_b), len_y, type_fpi_y, nbr_left, 5, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, 1, bt_b), len_y, type_fpi_y, nbr_left, 8, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! right nx+1(neibor 1) line inside effect zone, same tag as senders
                idx_global = i_start_global + nx + 1
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! receive fp_1,5,8 (nx+1) from right neighbor 
                    call MPI_Irecv(f_post(1, nx+1, bt_b), len_y, type_fpi_y, nbr_right, 1, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, nx+1, bt_b), len_y, type_fpi_y, nbr_right, 5, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, nx+1, bt_b), len_y, type_fpi_y, nbr_right, 8, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- 2 -(to left)-> nx+2 ---------------------*/
                ! left 2 line inside effect zone
                idx_global = i_start_global + 2
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! send fp_3,6,7 (2) to left neighbor
                    call MPI_Isend(f_post(3, 2, bt_b), len_y, type_fpi_y, nbr_left, 3, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(6, 2, bt_b), len_y, type_fpi_y, nbr_left, 6, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(7, 2, bt_b), len_y, type_fpi_y, nbr_left, 7, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! right nx+2(neibor 2) line inside effect zone, same tag as senders
                idx_global = i_start_global + nx + 2
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! receive fp_3,6,7 (nx+2) from right neighbor 
                    call MPI_Irecv(f_post(3, nx+2, bt_b), len_y, type_fpi_y, nbr_right, 3, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(6, nx+2, bt_b), len_y, type_fpi_y, nbr_right, 6, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(7, nx+2, bt_b), len_y, type_fpi_y, nbr_right, 7, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*----------message passing to right(i++)
                ! /*------------------- nx -(to right)-> 0 ---------------------*/
                ! right nx line inside effect zone
                idx_global = i_start_global + nx
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! send fp_3,6,7 (nx) to right neighbor
                    call MPI_Isend(f_post(3, nx, bt_b), len_y, type_fpi_y, nbr_right, 13, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(6, nx, bt_b), len_y, type_fpi_y, nbr_right, 16, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(7, nx, bt_b), len_y, type_fpi_y, nbr_right, 17, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! left 0 line inside effect zone
                idx_global = i_start_global
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! receive fp_3,6,7 (0) from left neighbor
                    call MPI_Irecv(f_post(3, 0, bt_b), len_y, type_fpi_y, nbr_left, 13, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(6, 0, bt_b), len_y, type_fpi_y, nbr_left, 16, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(7, 0, bt_b), len_y, type_fpi_y, nbr_left, 17, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- nx-1 -(to right)-> -1 ---------------------*/
                idx_global = i_start_global + nx - 1
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! send fp_1,5,8 (nx-1) to right neighbor
                    call MPI_Isend(f_post(1, nx-1, bt_b), len_y, type_fpi_y, nbr_right, 11, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, nx-1, bt_b), len_y, type_fpi_y, nbr_right, 15, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, nx-1, bt_b), len_y, type_fpi_y, nbr_right, 18, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! left -1 line inside effect zone
                idx_global = i_start_global - 1
                if (idx_global >= lf_b .AND. idx_global <= rg_b) then
                    ! receive fp_1,5,8 (-1) from left neighbor
                    call MPI_Irecv(f_post(1, -1, bt_b), len_y, type_fpi_y, nbr_left, 11, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, -1, bt_b), len_y, type_fpi_y, nbr_left, 15, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, -1, bt_b), len_y, type_fpi_y, nbr_left, 18, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif


                ! /*-------------------------- exchange message along y -----------------------------*/
                ! /*----------message passing to bottom(j--)
                ! /*------------------- 1 -(to bottom)-> ny+1 ---------------------*/
                ! bottom 1 line inside effect zone
                idx_global = j_start_global + 1
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! send fp_2,5,6 (1) to bottom neighbor
                    call MPI_Isend(f_post(2, lf_b, 1), len_x, type_fi_x, nbr_bottom, 22, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, lf_b, 1), len_x, type_fi_x, nbr_bottom, 25, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(6, lf_b, 1), len_x, type_fi_x, nbr_bottom, 26, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! top ny+1(neibor 1) line inside effect zone, same tag as senders
                idx_global = j_start_global + ny + 1
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! receive fp_2,5,6 (ny+1) lines from top neighbor 
                    call MPI_Irecv(f_post(2, lf_b, ny+1), len_x, type_fi_x, nbr_top, 22, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, lf_b, ny+1), len_x, type_fi_x, nbr_top, 25, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(6, lf_b, ny+1), len_x, type_fi_x, nbr_top, 26, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*------------------- 2 -(to bottom)-> ny+2 ---------------------*/
                ! bottom 2 line inside effect zone
                idx_global = j_start_global + 2
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! send fp_4,7,8 (2) to bottom neighbor
                    call MPI_Isend(f_post(4, lf_b, 2), len_x, type_fi_x, nbr_bottom, 24, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(7, lf_b, 2), len_x, type_fi_x, nbr_bottom, 27, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, lf_b, 2), len_x, type_fi_x, nbr_bottom, 28, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! top ny+2(neibor 2) line inside effect zone, same tag as senders
                idx_global = j_start_global + ny + 2
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! receive fp_4,7,8 (ny+1) lines from top neighbor 
                    call MPI_Irecv(f_post(4, lf_b, ny+2), len_x, type_fi_x, nbr_top, 24, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(7, lf_b, ny+2), len_x, type_fi_x, nbr_top, 27, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, lf_b, ny+2), len_x, type_fi_x, nbr_top, 28, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif


                ! /*----------message passing to top(j++)
                ! /*------------------- ny -(to top)-> 0 ---------------------*/
                ! top ny line inside effect zone
                idx_global = j_start_global + ny
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! send fp_4,7,8 (ny) to top neighbor
                    call MPI_Isend(f_post(4, lf_b, ny), len_x, type_fi_x, nbr_top, 34, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(7, lf_b, ny), len_x, type_fi_x, nbr_top, 37, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(8, lf_b, ny), len_x, type_fi_x, nbr_top, 38, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! bottom 0 line inside effect zone
                idx_global = j_start_global
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! receive fp_4,7,8 (0) from bottom neighbor
                    call MPI_Irecv(f_post(4, lf_b, 0), len_x, type_fi_x, nbr_bottom, 34, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(7, lf_b, 0), len_x, type_fi_x, nbr_bottom, 37, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(8, lf_b, 0), len_x, type_fi_x, nbr_bottom, 38, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                ! /*------------------- ny-1 -(to top)-> -1 ---------------------*/
                ! top ny-1 line inside effect zone
                idx_global = j_start_global + ny - 1
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! send fp_2,5,6 (ny-1) to top neighbor
                    call MPI_Isend(f_post(2, lf_b, ny-1), len_x, type_fi_x, nbr_top, 32, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(5, lf_b, ny-1), len_x, type_fi_x, nbr_top, 35, comm2d, req(cnt+2), rc)
                    call MPI_Isend(f_post(6, lf_b, ny-1), len_x, type_fi_x, nbr_top, 36, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif
                
                ! bottom -1 line inside effect zone
                idx_global = j_start_global
                if (idx_global >= bt_b .AND. idx_global <= tp_b) then
                    ! receive fp_2,5,6 (-1) from bottom neighbor
                    call MPI_Irecv(f_post(2, lf_b, -1), len_x, type_fi_x, nbr_bottom, 32, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(5, lf_b, -1), len_x, type_fi_x, nbr_bottom, 35, comm2d, req(cnt+2), rc)
                    call MPI_Irecv(f_post(6, lf_b, -1), len_x, type_fi_x, nbr_bottom, 36, comm2d, req(cnt+3), rc)
                    cnt = cnt + 3
                endif

                ! /*-------------------------- corner messages -----------------------------*/

        endif
    enddo

    call MPI_Waitall(cnt, req, MPI_STATUSES_IGNORE, rc)

end subroutine message_particle_f_post