    ! /* exchange messages along i direction*/
    ! // odd index 
    if (mod(block_x, 2) == 1) then  
        ! // messages passing along i++  --- send then receive
        if (block_x < nx_block) then    
            ! collect data to send (i=nx)
            do k = 1, nz
                do j = 1, ny
                    fps_send_i(1,j,k) = f_post(1,nx,j,k)
                    fps_send_i(2,j,k) = f_post(7,nx,j,k)
                    fps_send_i(3,j,k) = f_post(9,nx,j,k)
                    fps_send_i(4,j,k) = f_post(11,nx,j,k)
                    fps_send_i(5,j,k) = f_post(13,nx,j,k)
                enddo
            enddo
            call MPI_Ssend(fps_send_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank+1, 2, MPI_COMM_WORLD, rc)
        endif
        ! // receive messages along i++
        if (block_x < 1) then
            call MPI_Recv(fps_recv_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

            ! scatter the received data (i=0)
            do k = 1, nz 
                do j = 1, ny
                    f_post(1,0,j,k) = fps_recv_i(1,j,k)
                    f_post(7,0,j,k) = fps_recv_i(2,j,k)
                    f_post(9,0,j,k) = fps_recv_i(3,j,k)
                    f_post(11,0,j,k) = fps_recv_i(4,j,k)
                    f_post(13,0,j,k) = fps_recv_i(5,j,k)
                enddo
            enddo
        endif

        ! // messages passing along i--

    endif
    ! // even index 
    if (MOD(block_x, 2) == 0) then
        ! // messages passing along i++  --- receive then send
        ! // receive messages along i++
        if (block_x < 1) then
            call MPI_Recv(fps_recv_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

            ! scatter the received data (i=0)
            do k = 1, nz 
                do j = 1, ny
                    f_post(1,0,j,k) = fps_recv_i(1,j,k)
                    f_post(7,0,j,k) = fps_recv_i(2,j,k)
                    f_post(9,0,j,k) = fps_recv_i(3,j,k)
                    f_post(11,0,j,k) = fps_recv_i(4,j,k)
                    f_post(13,0,j,k) = fps_recv_i(5,j,k)
                enddo
            enddo
        endif

        ! // send messages along i++
        if (block_x < nx_block) then    
            ! collect data to send (i=nx)
            do k = 1, nz
                do j = 1, ny
                    fps_send_i(1,j,k) = f_post(1,nx,j,k)
                    fps_send_i(2,j,k) = f_post(7,nx,j,k)
                    fps_send_i(3,j,k) = f_post(9,nx,j,k)
                    fps_send_i(4,j,k) = f_post(11,nx,j,k)
                    fps_send_i(5,j,k) = f_post(13,nx,j,k)
                enddo
            enddo
            call MPI_Ssend(fps_send_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank+1, 2, MPI_COMM_WORLD, rc)
        endif

    endif