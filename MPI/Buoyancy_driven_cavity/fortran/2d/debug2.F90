subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    do j = 1, ny
        do i = 1, nx
            
    n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(1) = g(1,i,j)-g(3,i,j)
    n(2) = g(2,i,j)-g(4,i,j)
    n(3) = -4.0d0*g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)
        
            neq(0) = T(i,j)
            neq(1) = T(i,j)*u(i,j)
            neq(2) = T(i,j)*v(i,j)
            neq(3) = T(i,j)*paraA
            neq(4) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qnu
            q(4) = Qnu
        
            do alpha=0,4
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j) = 0.2d0*n_post(0)-0.2d0*n_post(3)
    g_post(1,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(2,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
    g_post(3,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(4,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collisionT


subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    do j = 1, ny
        do i = 1, nx
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                
                    g(alpha,i,j) = g_post(alpha,ip,jp)
                enddo
        enddo
    enddo
    
    return
end subroutine streamingT


subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0,y0,q,Cd1,Cd2,Cd3,Cd4


#ifdef HorizontalWallsAdiabatic
    !Bottom side
    if (coords(1) == 0) then
        do i = 1, nx 
            g(2, i, 1) = g_post(4, i, 1)
        enddo
    endif

    ! Top side
    if (coords(1) == dims(1) - 1) then
        do i = 1, nx 
            g(4, i, ny) = g_post(2, i, ny)
        enddo
    endif
#endif


#ifdef HorizontalWallsConstT
   !Bottom side
    if (coords(1) == 0) then
        do i = 1, nx 
            g(2, i, 1) = -g_post(4, i, 1)+(4.0d0+paraA)/10.0d0*Thot
        enddo
    endif

    ! Top side
    if (coords(1) == dims(1) - 1) then
        do i = 1, nx 
            g(4, i, ny) = -g_post(2, i, ny)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
    endif
#endif

#ifdef VerticalWallsAdiabatic
    !Left side
    if (coords(0) == 0) then
        do j = 1, ny 
            g(1,1,j) = g_post(3,1,j)
        enddo
    endif

    !Right side
    if (coords(0) == dims(0) - 1) then
        do j = 1, ny 
            g(3,nx,j) = g_post(1,nx,j)
        enddo
    endif
#endif

#ifdef VerticalWallsConstT
    !Left side
    if (coords(0) == 0) then
        do j = 1, ny 
            g(1,1,j) = -g_post(3,1,j)+(4.0d0+paraA)/10.0d0*Thot
        enddo
    endif

    !Right side
    if (coords(0) == dims(0) - 1) then
        do j = 1, ny 
            g(3,nx,j) = -g_post(1,nx,j)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
    endif
#endif


    
! #ifdef VerticalWallsPeriodicalT
!     !$omp parallel do default(none) shared(g,g_post) private(j) 
!     do j=1,ny
!         !Left side
!         g(1,1,j) = g_post(1,nx,j)

!         !Right side
!         g(3,nx,j) = g_post(3,1,j)
!     enddo
!     !$omp end parallel do
! #endif

    return
end subroutine bouncebackT



subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: a
    
    do j = 1, ny
        do i = 1, nx
                T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        enddo
    enddo

    return
end subroutine macroT
    

#ifdef steadyFlow
    subroutine check()
    use mpi
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    

    do j=1,ny
        do i=1,nx
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                error5 = error5+dABS( T(i,j)-Tp(i,j) )
                error6 = error6+dABS( T(i,j) )
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j)
                Tp(i,j) = T(i,j)
        enddo
    enddo

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        open(unit=01,file='convergence.log',status='unknown',position='append')
            write(01,*) itc,' ',errorU,' ',errorT
        close(01)
        
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    
    return
    end subroutine check
#endif



subroutine message_passing_f()
    use mpi
    use commondata
    implicit none

    ! message tag:  0, 1, 2, discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Sendrecv(f_post(2, 1, ny), 1, f_row_x, nbr_top, 2, &
                    f_post(2, 1, 0), 1, f_row_x, nbr_bottom, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, 1, ny), 1, f_row_x, nbr_top, 5, &
                    f_post(5, 1, 0), 1, f_row_x, nbr_bottom, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(6, 1, ny), 1, f_row_x, nbr_top, 6, &
                    f_post(6, 1, 0), 1, f_row_x, nbr_bottom, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to bottom(j--)
    call MPI_Sendrecv(f_post(4, 1, 1), 1, f_row_x, nbr_bottom, 4, &
                    f_post(4, 1, ny+1), 1, f_row_x, nbr_top, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(7, 1, 1), 1, f_row_x, nbr_bottom, 7, &
                    f_post(7, 1, ny+1), 1, f_row_x, nbr_top, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, 1, 1), 1, f_row_x, nbr_bottom, 8, &
                    f_post(8, 1, ny+1), 1, f_row_x, nbr_top, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f_post(1, nx, 1), 1, f_column_y, nbr_right, 1, &
                    f_post(1, 0, 1), 1, f_column_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, nx, 1), 1, f_column_y, nbr_right, 5, &
                    f_post(5, 0, 1), 1, f_column_y, nbr_left, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, nx, 1), 1, f_column_y, nbr_right, 8, &
                    f_post(8, 0, 1), 1, f_column_y, nbr_left, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    ! message passing to left(i--)
    call MPI_Sendrecv(f_post(3, 1, 1), 1, f_column_y, nbr_left, 3, &
                    f_post(3, nx+1, 1), 1, f_column_y, nbr_right, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)       

    call MPI_Sendrecv(f_post(6, 1, 1), 1, f_column_y, nbr_left, 6, &
                    f_post(6, nx+1, 1), 1, f_column_y, nbr_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    call MPI_Sendrecv(f_post(7, 1, 1), 1, f_column_y, nbr_left, 7, &
                    f_post(7, nx+1, 1), 1, f_column_y, nbr_right, 7, &
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
                

end subroutine message_passing_f





subroutine message_passing_g()
    use mpi
    use commondata
    implicit none

    ! message tag:  0, 1, 2, discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Sendrecv(g_post(2, 1, ny), 1, g_row_x, nbr_top, 2, &
                    g_post(2, 1, 0), 1, g_row_x, nbr_bottom, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)


    ! message passing to bottom(j--)
    call MPI_Sendrecv(g_post(4, 1, 1), 1, g_row_x, nbr_bottom, 4, &
                    g_post(4, 1, ny+1), 1, g_row_x, nbr_top, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)


    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(g_post(1, nx, 1), 1, g_column_y, nbr_right, 1, &
                    g_post(1, 0, 1), 1, g_column_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to left(i--)
    call MPI_Sendrecv(g_post(3, 1, 1), 1, g_column_y, nbr_left, 3, &
                    g_post(3, nx+1, 1), 1, g_column_y, nbr_right, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)       
                

end subroutine message_passing_g



subroutine output()
    use mpi
    use commondata
    integer :: i, j
    integer :: p_rank, num(0:3) ,dx = 0, dy = 0, new_coords(0:1)
    real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :), total_T(:, :), stream(:, :), vorticity(:, :)
    real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :), tmp_T(:, :)

    if (rank2d > 0) then
        ! rank != 0 send data
        num(0) = nx
        num(1) = ny
        num(2) = i_start_global
        num(3) = j_start_global
        ! send to rank 0
        call MPI_Send(num, 4, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size and origion
        call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
        call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
        call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
        call MPI_Send(T, nx*ny, MPI_REAL8, 0, 4, comm2d, rc)
    else
        ! rank 0 collect data
        allocate(total_u(total_nx, total_ny))
        allocate(total_v(total_nx, total_ny))
        allocate(total_rho(total_nx, total_ny))
        allocate(total_T(total_nx, total_ny))

        dx = i_start_global
        dy = j_start_global

        ! collect data from rank 0
        do j = 1, ny
            do i = 1, nx
                total_u(dx + i, dy + j) = u(i, j)
                total_v(dx + i, dy + j) = v(i, j)
                total_rho(dx + i, dy + j) = rho(i, j)
                total_T(dx + i, dy + j) = T(i, j)
            enddo
        enddo

         ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) - 1
            ! receive the block size and origion
            call MPI_Recv(num, 4, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)
            ! creat buffer
            allocate(tmp_u(num(0), num(1)))
            allocate(tmp_v(num(0), num(1)))
            allocate(tmp_rho(num(0), num(1)))
            allocate(tmp_T(num(0), num(1)))
            
            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_T, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 4, comm2d, MPI_STATUS_IGNORE, rc)

            ! determine the origin
            dx = num(2)
            dy = num(3)

            ! assign data
            do j = 1, num(1)
                do i = 1, num(0)
                    total_u(dx + i, dy + j) = tmp_u(i, j)
                    total_v(dx + i, dy + j) = tmp_v(i, j)
                    total_rho(dx + i, dy + j) = tmp_rho(i, j)
                    total_T(dx + i, dy + j) = tmp_T(i, j)
                enddo
            enddo

            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_rho)
            deallocate(tmp_T)
        enddo


        allocate(stream(total_nx, total_ny))
        allocate(vorticity(total_nx, total_ny))

        call compute_stream_vorticity(stream, vorticity, total_u, total_v, total_nx, total_ny)
        call output_Tecplot(xp, yp, total_u, total_v, total_rho, total_T, stream, vorticity, total_nx, total_ny, itc)
        call output_binary(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        call output_ASCII(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        call out_Velocity_Nu(total_u, total_v, total_T, total_nx, total_ny, diffusivity, lengthUnit)

        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_rho)
        deallocate(total_T)
        deallocate(stream)
        deallocate(vorticity)

    endif


end subroutine output

subroutine out_Velocity_Nu(u, v, T, nx, ny, diffusivity, lengthUnit)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: diffusivity, lengthUnit
    real(8), intent(in) :: u(nx, ny), v(nx, ny), T(nx, ny)
    integer :: i, j, nx_half, ny_half
    real(8) :: max_u, max_v, u_i, v_i
    real(8) :: Nu, Nu_vol, Nu_0, Nu_half, Nu_max, Nu_min, Nu_max_id, Nu_min_id

    nx_half = (1 + nx) / 2
    ny_half = (1 + ny) / 2

    ! find maximun u along vertical mid-plane, and the maximun v along horizontal mid-plane
    
    max_v = -1e30
    do i = 1, nx
        if (v(i, ny_half) > max_v) then
            max_v = v(i, ny_half)
            v_i = i
        endif
    enddo

    max_u = -1e30
    do j = 1, ny
        if (u(nx_half, j) > max_u) then
            max_u = u(nx_half, j)
            u_i = j
        endif
    enddo

    max_u = max_u / (diffusivity / lengthUnit)
    max_v = max_v / (diffusivity / lengthUnit)
    u_i = u_i / ny
    v_i = v_i / nx

    open(unit=77,file="result_parameters.txt",status='unknown',position="append")
        write(77,300) max_u, u_i, max_v, v_i
        write(*,300) max_u, u_i, max_v, v_i
300     format('max_u = ', f13.7," at y = ", f13.7, '  max_v = ', f13.7, " at x = ", f13.7)
    close(77)

    

    ! calculate Nu 
    Nu_vol = 0.0d0
    Nu_0 = 0.0d0
    Nu_half = 0.0d0
    Nu_max = -1e30
    Nu_min = 1e30
    do j = 1, ny
        do i = 1, nx
            Nu = v(i, j) * T(i, j)
            Nu_vol = Nu_vol + Nu
            if (i == 1) then     !! hot wall
                Nu_0 = Nu_0 + Nu
                if (Nu > Nu_max) then
                    Nu_max = Nu
                    Nu_max_id = j
                endif
                if (Nu < Nu_min) then
                    Nu_min = Nu
                    Nu_min_id = j
                endif
            endif

            if (i == nx_half) then  !! vertical mid-plane
                Nu_half = Nu_half + Nu
            endif
        enddo
    enddo

    Nu_vol = Nu_vol / (nx * ny) * lengthUnit / diffusivity + 1.0
    Nu_0 = Nu_0 / ny * lengthUnit / diffusivity + 1.0
    Nu_half = Nu_half / ny * lengthUnit / diffusivity + 1.0
    Nu_max = Nu_max * lengthUnit / diffusivity + 1.0
    Nu_min = Nu_min * lengthUnit / diffusivity + 1.0
    Nu_max_id = Nu_max_id / ny
    Nu_min_id = Nu_min_id / ny


    open(unit=77,file="result_parameters.txt",status='unknown',position="append")
        write(77,301) Nu_vol, Nu_0, Nu_half
        write(77,302) Nu_max, Nu_max_id, Nu_min
        write(*,301) Nu_vol, Nu_0, Nu_half
        write(*,302) Nu_max, Nu_max_id, Nu_min, Nu_min_id
301     format('Nu_vol = ', f13.7, '  Nu_0 = ', f13.7, '  Nu_half = ', f13.7)
302     format('Nu_max = ', f13.7," at y = ", f13.7, '  Nu_min = ', f13.7, " at y = ", f13.7)
    close(77)


end subroutine out_Velocity_Nu


subroutine output_binary(u, v, rho, T, nx, ny, itc)
    use ioFolder
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny)
    integer(kind=4) :: i, j
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc
#endif
#ifdef unsteadyFlow
    binFileNum = binFileNum+1
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif
    filename = adjustl(filename)

    open(unit=03,file=trim(binFolderPrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(03) ((u(i,j),i=1,nx),j=1,ny)
    write(03) ((v(i,j),i=1,nx),j=1,ny)
    write(03) ((T(i,j),i=1,nx),j=1,ny)
    close(03)

    return
end subroutine output_binary



subroutine output_Tecplot(xp, yp, u, v, rho, T, stream, vorticity, nx, ny, itc)
    use ioFolder
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny), stream(nx, ny), vorticity(nx, ny)
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7,V8
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1
    write(filename,'(i12.12)') pltFileNum
#endif
    filename = adjustl(filename)
    
    open(unit=41,file=trim(pltFolderPrefix)//"-"//trim(filename)//'.plt',form='Unformatted')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--Integer value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file

    write(41) 8

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)
    V6='rho'
    call dumpstring(V6)
    V7='stream'
    call dumpstring(V7)
    V8='vorticity'
    call dumpstring(V8)


    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName='ZONE 001'
    call dumpstring(zoneName)

    !---------Zone Color
    write(41) -1

    !---------ZoneType
    write(41) 0

    !---------DataPacking 0=Block, 1=Point
    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax
    write(41) nx
    write(41) ny
    write(41) kmax

    !-----------1=Auxiliary name/value pair to follow
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    !----zone ------------------------------------------------------------
    write(41) zoneMarker

    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(u(i,j))
                write(41) real(v(i,j))
                write(41) real(T(i,j))
                write(41) real(rho(i,j))
                write(41) real(stream(i,j))
                write(41) real(vorticity(i,j))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
end subroutine output_Tecplot


subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
end subroutine dumpstring



! subroutine backupData()
!     use commondata
!     implicit none
!     integer(kind=4) :: i, j, alpha
!     character(len=100) :: filename

! #ifdef steadyFlow
!     write(filename,*) itc
! #endif
! #ifdef unsteadyFlow
!     if(loadInitField.EQ.0) write(filename,*) binFileNum
!     if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
! #endif
!     filename = adjustl(filename)

!     open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
!     write(05) (((f(alpha,i,j),alpha=0,8), i=1,nx), j=1,ny)
!     write(05) (((g(alpha,i,j),alpha=0,4), i=1,nx), j=1,ny)
!     write(05) ((u(i,j), i=1,nx), j=1,ny)
!     write(05) ((v(i,j), i=1,nx), j=1,ny)
!     write(05) ((T(i,j), i=1,nx), j=1,ny)
!     close(05)
    
!     open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
!     write(00,*) "Backup  f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
!     close(00)

!     return
! end subroutine backupData


subroutine compute_stream_vorticity(stream, vorticity, u, v, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: u(nx, ny), v(nx, ny)
    real(8), intent(out) :: stream(nx, ny), vorticity(nx, ny)
    integer :: i, j


    do j = 1, ny
        stream(1, j) = -v(1, j) / 4.0d0
        stream(2, j) = -3.0d0 / 4.0d0 * v(1, j) - v(2, j) / 2.0d0
        do i = 2, nx-1
            stream(i+1, j) = stream(i-1, j) - (v(i-1, j) + 4.0d0 * v(i, j) + v(i+1, j)) / 3.0d0
        enddo
    enddo

    do j = 2, ny-1
        do i = 2, nx-1
            vorticity(i, j) = 1.0d0 / 3.0d0 * (v(i+1, j) - v(i-1, j)) &
                                + 1.0d0 / 12.0d0 * (v(i+1, j+1) - v(i-1, j-1)) &
                                + 1.0d0 / 12.0d0 * (v(i+1, j-1) - v(i-1, j+1)) &
                            - 1.0d0 / 3.0d0 * (u(i, j+1) - u(i, j-1)) &
                                - 1.0d0 / 12.0d0 * (u(i+1, j+1) - u(i-1, j-1)) &
                                - 1.0d0 / 12.0d0 * (u(i-1, j+1) - u(i+1, j-1))
        enddo
    enddo

    do j = 2, ny-1
        vorticity(1, j) = 2.0d0 * v(2, j)
        vorticity(nx, j) = -2.0d0 * v(nx-1, j)
    enddo

    do i = 2, nx-1
        vorticity(i, 1) = 2.0d0 * u(i, 2)
        vorticity(i, ny) = -2.0d0 * u(i, ny-1)
    enddo

    vorticity(1, ny) = 0
    vorticity(nx, ny) = 0

    vorticity(1, 1) = 0
    vorticity(nx, 1) = 0

end subroutine compute_stream_vorticity



subroutine output_ASCII(u, v, rho, T, nx, ny, itc)
    use ioFolder
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny)
    integer(kind=4) :: i, j
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc
#endif
#ifdef unsteadyFlow
    binFileNum = binFileNum+1
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif
    filename = adjustl(filename)

    open(unit=73,file="data"//"-"//trim(filename)//'.dat', status='unknown')
        write(73, 730) ((u(i,j),i=1,nx),j=1,ny)
        write(73, 730) ((v(i,j),i=1,nx),j=1,ny)
        write(73, 730) ((T(i,j),i=1,nx),j=1,ny)
    close(73)

730     format((e13.6,' '))
    return
end subroutine output_ASCII