subroutine message_passing()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    ! surface data to send:   fps_send_i(5,j,k) is data at surface i=1 and i = nx
    real(8) :: fps_send_i(5,ny,nz), fps_send_j(5,nx,nz), fps_send_k(5,nx,ny)
    ! surface data to receive
    real(8) :: fps_recv_i(5,ny,nz), fps_recv_j(5,nx,nz), fps_recv_k(5,nx,ny)

    ! line data to send  ---  fpl_send_z(k) is the line parallel to axis z
    real(8) :: fpl_send_z(nz), fpl_send_y(ny), fpl_send_x(nx)
    real(8) :: fpl_recv_z(nz), fpl_recv_y(ny), fpl_recv_x(nx)

    ! we don't need to exchange fp for cornor points


    !!! /* ------------- exchange message with surfaces ----------------- 
    !!!    message tag(send direction):  1 -> i--  /  2 -> i++   /  3 -> j--  /
    !!!                                  4 -> j++  /  5 -> k--   /  6 -> k=++  */


    ! /* exchange messages along i direction*/
    ! // odd index 
    if (mod(block_x, 2) == 1) then  
        ! /* messages passing along i-- ---- send then receive */
        if (block_x > 1) then
            call MPI_Send_Surface(f_post, rank-1, 1, nx, ny, nz)
        endif
        if (block_x < nx_block) then
            call MPI_Recv_Surface(f_post, rank+1, 1, nx, ny, nz)
        endif

        ! /* messages passing along i++  --- send then receive */
        if (block_x < nx_block) then    
            call MPI_Send_Surface(f_post, rank+1, 2, nx, ny, nz)
        endif
        if (block_x > 1) then
            call MPI_Recv_Surface(f_post, rank-1, 2, nx, ny, nz)
        endif
    endif

    ! // even index 
    if (MOD(block_x, 2) == 0) then
        ! /* messages passing along i-- ---- receive then send */
        if (block_x < nx_block) then
            call MPI_Recv_Surface(f_post, rank+1, 1, nx, ny, nz)
        endif
        if (block_x > 1) then
            call MPI_Send_Surface(f_post, rank-1, 1, nx, ny, nz)
        endif

        ! // messages passing along i++  --- receive then send
        if (block_x > 1) then
            call MPI_Recv_Surface(f_post, rank-1, 2, nx, ny, nz)
        endif
        if (block_x < nx_block) then 
            call MPI_Send_Surface(f_post, rank+1, 2, nx, ny, nz)   
        endif
    endif












    ! /* exchange messages along k direction*/
    ! // odd index 
    if (mod(block_z, 2) == 1) then  
        ! /* messages passing along k++  --- send then receive */
        ! // send messages along k++
        if (block_z < nz_block) then    
            ! ! collect data to send (k=nz)
            do j = 1, ny
                do i = 1, nx
                    fps_send_k(1,i,j) = f_post(5,i,j,nz)
                    fps_send_k(2,i,j) = f_post(11,i,j,nz)
                    fps_send_k(3,i,j) = f_post(12,i,j,nz)
                    fps_send_k(4,i,j) = f_post(15,i,j,nz)
                    fps_send_k(5,i,j) = f_post(16,i,j,nz)
                enddo
            enddo
            call MPI_Send(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block, 6, MPI_COMM_WORLD, rc)
        endif

        ! // receive messages along k++
        if (block_z > 1) then
            call MPI_Recv(fps_recv_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
            ! scatter the received data (k=0)
            do j = 1, ny
                do i = 1, nx
                    f_post(5,i,j,0) = fps_recv_k(1,i,j)
                    f_post(11,i,j,0) = fps_recv_k(2,i,j)
                    f_post(12,i,j,0) = fps_recv_k(3,i,j)
                    f_post(15,i,j,0) = fps_recv_k(4,i,j)
                    f_post(16,i,j,0) = fps_recv_k(5,i,j)
                enddo
            enddo
        endif

        ! /* messages passing along k-- ---- send then receive */
        ! // send messages along k--
        if (block_z > 1) then
            ! collect data to send (k=1)
            do j = 1, ny
                do i = 1, nx
                    fps_send_k(1,i,j) = f_post(6,i,j,1)
                    fps_send_k(2,i,j) = f_post(13,i,j,1)
                    fps_send_k(3,i,j) = f_post(14,i,j,1)
                    fps_send_k(4,i,j) = f_post(17,i,j,1)
                    fps_send_k(5,i,j) = f_post(18,i,j,1)
                enddo
            enddo
            call MPI_Send(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block, 5, MPI_COMM_WORLD, rc)
        endif

        ! // receive messages along k--
        if (block_z < nz_block) then
            call MPI_Recv(fps_recv_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)                
            ! scatter the received data(k=nz+1)
            do j = 1, ny
                do i = 1, nx
                    f_post(6,i,j,nz+1) = fps_recv_k(1,i,j)
                    f_post(13,i,j,nz+1) = fps_recv_k(2,i,j)
                    f_post(14,i,j,nz+1) = fps_recv_k(3,i,j)
                    f_post(17,i,j,nz+1) = fps_recv_k(4,i,j)
                    f_post(18,i,j,nz+1) = fps_recv_k(5,i,j)
                enddo
            enddo
        endif
    endif

    ! // even index 
    if (mod(block_z, 2) == 0) then  
        ! /* messages passing along k++  --- receive then send */
        ! // receive messages along k++
        if (block_z > 1) then
            call MPI_Recv(fps_recv_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)                
            ! scatter the received data (k=0)
            do j = 1, ny
                do i = 1, nx
                    f_post(5,i,j,0) = fps_recv_k(1,i,j)
                    f_post(11,i,j,0) = fps_recv_k(2,i,j)
                    f_post(12,i,j,0) = fps_recv_k(3,i,j)
                    f_post(15,i,j,0) = fps_recv_k(4,i,j)
                    f_post(16,i,j,0) = fps_recv_k(5,i,j)
                enddo
            enddo
        endif

        ! // send messages along k++
        if (block_z < nz_block) then    
            ! collect data to send (k=nz)
            do j = 1, ny
                do i = 1, nx
                    fps_send_k(1,i,j) = f_post(5,i,j,nz)
                    fps_send_k(2,i,j) = f_post(11,i,j,nz)
                    fps_send_k(3,i,j) = f_post(12,i,j,nz)
                    fps_send_k(4,i,j) = f_post(15,i,j,nz)
                    fps_send_k(5,i,j) = f_post(16,i,j,nz)
                enddo
            enddo
            call MPI_Send(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block, 6, MPI_COMM_WORLD, rc)
        endif

        ! /* messages passing along k-- ---- receive then send */
        ! // receive messages along k--
        if (block_z < nz_block) then
            call MPI_Recv(fps_recv_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)         
            ! scatter the received data(k=nz+1)
            do j = 1, ny
                do i = 1, nx
                    f_post(6,i,j,nz+1) = fps_recv_k(1,i,j)
                    f_post(13,i,j,nz+1) = fps_recv_k(2,i,j)
                    f_post(14,i,j,nz+1) = fps_recv_k(3,i,j)
                    f_post(17,i,j,nz+1) = fps_recv_k(4,i,j)
                    f_post(18,i,j,nz+1) = fps_recv_k(5,i,j)
                enddo
            enddo
        endif
        ! // send messages along k--
        if (block_z > 1) then
            ! collect data to send (k=1)
            do j = 1, ny
                do i = 1, nx
                    fps_send_k(1,i,j) = f_post(6,i,j,1)
                    fps_send_k(2,i,j) = f_post(13,i,j,1)
                    fps_send_k(3,i,j) = f_post(14,i,j,1)
                    fps_send_k(4,i,j) = f_post(17,i,j,1)
                    fps_send_k(5,i,j) = f_post(18,i,j,1)
                enddo
            enddo
            call MPI_Send(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block, 5, MPI_COMM_WORLD, rc)
        endif
    endif



end subroutine message_passing



subroutine MPI_Send_Surface(f_post, tar_rank, tag, nx, ny, nz)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, nz, tar_rank, tag
    real(8), intent(in) :: f_post(0:18, 0:nx+1, 0:ny+1, 0:nz+1)
    integer :: i, j, k, rc
    real(8),allocatable :: buffer(:, :, :)

    !!! /* ------------- exchange message with surfaces ----------------- 
    !!!    message tag(send direction):  1 -> i--  /  2 -> i++   /  3 -> j--  /
    !!!                                  4 -> j++  /  5 -> k--   /  6 -> k=++  */\

    message_tag: select case (tag)
    case (1)
        allocate(buffer(5, ny, nz))
        ! send to i--, collect data at (i=1)
        do k = 1, nz
            do j = 1, ny
                buffer(1,j,k) = f_post(2,1,j,k)
                buffer(2,j,k) = f_post(8,1,j,k)
                buffer(3,j,k) = f_post(10,1,j,k)
                buffer(4,j,k) = f_post(12,1,j,k)
                buffer(5,j,k) = f_post(14,1,j,k)
            enddo
        enddo
        call MPI_Send(buffer, 5*ny*nz, MPI_DOUBLE_PRECISION, tar_rank, 1, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case (2)
        allocate(buffer(5, ny, nz))
        ! send to i++, collect data at (i=nx)
        do k = 1, nz
            do j = 1, ny
                buffer(1,j,k) = f_post(1,nx,j,k)
                buffer(2,j,k) = f_post(7,nx,j,k)
                buffer(3,j,k) = f_post(9,nx,j,k)
                buffer(4,j,k) = f_post(11,nx,j,k)
                buffer(5,j,k) = f_post(13,nx,j,k)
            enddo
        enddo
        call MPI_Send(buffer, 5*ny*nz, MPI_DOUBLE_PRECISION, tar_rank, 2, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case (3)
        allocate(buffer(5, nx, nz))
        ! send to j--, collect data at (j=1)
        do k = 1, nz
            do i = 1, nx
                buffer(1,i,k) = f_post(4,i,1,k)
                buffer(2,i,k) = f_post(9,i,1,k)
                buffer(3,i,k) = f_post(10,i,1,k)
                buffer(4,i,k) = f_post(16,i,1,k)
                buffer(5,i,k) = f_post(18,i,1,k)
            enddo
        enddo
        call MPI_Send(buffer, 5*nx*nz, MPI_DOUBLE_PRECISION, tar_rank, 3, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case (4)
        allocate(buffer(5, nx, nz))
        ! send to j++, collect data at (j=1)
        do k = 1, nz
            do i = 1, nx
                buffer(1,i,k) = f_post(3,i,ny,k)
                buffer(2,i,k) = f_post(7,i,ny,k)
                buffer(3,i,k) = f_post(8,i,ny,k)
                buffer(4,i,k) = f_post(15,i,ny,k)
                buffer(5,i,k) = f_post(17,i,ny,k)
            enddo
        enddo
        call MPI_Send(buffer, 5*nx*nz, MPI_DOUBLE_PRECISION, tar_rank, 4, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case (5)
        allocate(buffer(5, nx, ny))
        ! send to k--, collect data to send (k=1)
        do j = 1, ny
            do i = 1, nx
                buffer(1,i,j) = f_post(6,i,j,1)
                buffer(2,i,j) = f_post(13,i,j,1)
                buffer(3,i,j) = f_post(14,i,j,1)
                buffer(4,i,j) = f_post(17,i,j,1)
                buffer(5,i,j) = f_post(18,i,j,1)
            enddo
        enddo
        call MPI_Send(buffer, 5*nx*ny, MPI_DOUBLE_PRECISION, tar_rank, 5, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case (6)
        allocate(buffer(5, nx, ny))
        ! send to k++, collect data to send (k=nz)
        do j = 1, ny
            do i = 1, nx
                buffer(1,i,j) = f_post(5,i,j,nz)
                buffer(2,i,j) = f_post(11,i,j,nz)
                buffer(3,i,j) = f_post(12,i,j,nz)
                buffer(4,i,j) = f_post(15,i,j,nz)
                buffer(5,i,j) = f_post(16,i,j,nz)
            enddo
        enddo
        call MPI_Send(buffer, 5*nx*ny, MPI_DOUBLE_PRECISION, tar_rank, 6, MPI_COMM_WORLD, rc)
        deallocate(buffer)
    case default
        write(*, *) "Send surface data -- wrong message tag, make sure 1 <= tag <= 6."
    end select message_tag

end subroutine MPI_Send_Surface



subroutine MPI_Recv_Surface(f_post, tar_rank, tag, nx, ny, nz)
    use mpi
    implicit none
    integer, intent(in) :: nx, ny, nz, tar_rank, tag
    real(8), intent(inout) :: f_post(0:18, 0:nx+1, 0:ny+1, 0:nz+1)
    integer :: i, j, k, rc
    real(8),allocatable :: buffer(:, :, :)

    !!! /* ------------- exchange message with surfaces ----------------- 
    !!!    message tag(send direction):  1 -> i--  /  2 -> i++   /  3 -> j--  /
    !!!                                  4 -> j++  /  5 -> k--   /  6 -> k=++  */\

    message_tag: select case (tag)
    case (1)
        allocate(buffer(5, ny, nz))
        ! receive from i--
        call MPI_Recv(buffer, 5*ny*nz, MPI_DOUBLE_PRECISION, tar_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (i=nx+1)
        do k = 1, nz
            do j = 1, ny
                f_post(2,nx+1,j,k) = buffer(1,j,k)
                f_post(8,nx+1,j,k) = buffer(2,j,k)
                f_post(10,nx+1,j,k) = buffer(3,j,k)
                f_post(12,nx+1,j,k) = buffer(4,j,k)
                f_post(14,nx+1,j,k) = buffer(5,j,k)
            enddo
        enddo
        deallocate(buffer)
    case (2)
        allocate(buffer(5, ny, nz))
        ! receive from i++
        call MPI_Recv(buffer, 5*ny*nz, MPI_DOUBLE_PRECISION, tar_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (i=0)
        do k = 1, nz 
            do j = 1, ny
                f_post(1,0,j,k) = buffer(1,j,k)
                f_post(7,0,j,k) = buffer(2,j,k)
                f_post(9,0,j,k) = buffer(3,j,k)
                f_post(11,0,j,k) = buffer(4,j,k)
                f_post(13,0,j,k) = buffer(5,j,k)
            enddo
        enddo
        deallocate(buffer)
    case (3)
        allocate(buffer(5, nx, nz))
        ! receive from j--
        call MPI_Recv(buffer, 5*nx*nz, MPI_DOUBLE_PRECISION, tar_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (j=ny+1)
        do k = 1, nz
            do i = 1, nx
                f_post(4,i,ny+1,k) = buffer(1,i,k)
                f_post(9,i,ny+1,k) = buffer(2,i,k)
                f_post(10,i,ny+1,k) = buffer(3,i,k)
                f_post(16,i,ny+1,k) = buffer(4,i,k)
                f_post(18,i,ny+1,k) = buffer(5,i,k)
            enddo
        enddo
        deallocate(buffer)
    case (4)
        allocate(buffer(5, nx, nz))
        ! receive from j++
        call MPI_Recv(buffer, 5*nx*nz, MPI_DOUBLE_PRECISION, tar_rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (j=0)
        do k = 1, nz
            do i = 1, nx
                f_post(3,i,0,k) = buffer(1,i,k)
                f_post(7,i,0,k) = buffer(2,i,k)
                f_post(8,i,0,k) = buffer(3,i,k)
                f_post(15,i,0,k) = buffer(4,i,k)
                f_post(17,i,0,k) = buffer(5,i,k)
            enddo
        enddo
        deallocate(buffer)
    case (5)
        allocate(buffer(5, nx, ny))
        ! receive from k--
        call MPI_Recv(buffer, 5*nx*ny, MPI_DOUBLE_PRECISION, tar_rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (k=ny+1)
        do j = 1, ny
            do i = 1, nx
                f_post(6,i,j,nz+1) = buffer(1,i,j)
                f_post(13,i,j,nz+1) = buffer(2,i,j)
                f_post(14,i,j,nz+1) = buffer(3,i,j)
                f_post(17,i,j,nz+1) = buffer(4,i,j)
                f_post(18,i,j,nz+1) = buffer(5,i,j)
            enddo
        enddo
        deallocate(buffer)
    case (6)
        allocate(buffer(5, nx, ny))
        ! receive from k++
        call MPI_Recv(buffer, 5*nx*ny, MPI_DOUBLE_PRECISION, tar_rank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        ! scatter data at (k=0)
        do j = 1, ny
            do i = 1, nx
                f_post(5,i,j,0) = buffer(1,i,j)
                f_post(11,i,j,0) = buffer(2,i,j)
                f_post(12,i,j,0) = buffer(3,i,j)
                f_post(15,i,j,0) = buffer(4,i,j)
                f_post(16,i,j,0) = buffer(5,i,j)
            enddo
        enddo
        deallocate(buffer)
    case default
        write(*, *) "Send surface data -- wrong message tag, make sure 1 <= tag <= 6."
    end select message_tag

end subroutine MPI_Recv_Surface