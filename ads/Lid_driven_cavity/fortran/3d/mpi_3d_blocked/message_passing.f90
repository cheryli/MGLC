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

    integer :: buffer_size, max_size, s1
    character, allocatable :: buffer(:)

    max_size = 5 * max(max(nx*ny, nx*nz), ny*nz)

    call MPI_Pack_size(max_size, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, s1, rc)

    buffer_size = 2 * MPI_BSEND_OVERHEAD + s1;

    allocate(buffer(buffer_size))

    call MPI_Buffer_attach(buffer, buffer_size, rc)


    !!! /* ------------- exchange message with surfaces ----------------- 
    !!!    message tag(for sender):  1 -> i=1  /  2 -> i=nx   /  3 -> j=1  /
    !!!                              4 -> j=ny   /   5 -> k=1    /   6 -> k=nz  */

    ! /* exchange message with i=1 surface --- send then receive*/
    if (block_x > 1) then
        ! collect data to send (i=1)
        do k = 1, nz
            do j = 1, ny
                fps_send_i(1,j,k) = f_post(2,1,j,k)
                fps_send_i(2,j,k) = f_post(8,1,j,k)
                fps_send_i(3,j,k) = f_post(10,1,j,k)
                fps_send_i(4,j,k) = f_post(12,1,j,k)
                fps_send_i(5,j,k) = f_post(14,1,j,k)
            enddo
        enddo

        call MPI_Bsend(fps_send_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank-1, 1, MPI_COMM_WORLD, rc)
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

    ! /* exchange message with i=nx surface --- receive then send*/
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

        call MPI_Bsend(fps_send_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank+1, 2, MPI_COMM_WORLD, rc)
        call MPI_Recv(fps_recv_i, 5*ny*nz, MPI_DOUBLE_PRECISION, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        
        ! scatter the received data (i=nx+1)
        do k = 1, nz
            do j = 1, ny
                f_post(2,nx+1,j,k) = fps_recv_i(1,j,k)
                f_post(8,nx+1,j,k) = fps_recv_i(2,j,k)
                f_post(10,nx+1,j,k) = fps_recv_i(3,j,k)
                f_post(12,nx+1,j,k) = fps_recv_i(4,j,k)
                f_post(14,nx+1,j,k) = fps_recv_i(5,j,k)
            enddo
        enddo
    endif

    ! /* exchange message with j=1 surface --- send then receive*/
    if (block_y > 1) then
        ! collect data to send (j=1)
        do k = 1, nz
            do i = 1, nx
                fps_send_j(1,i,k) = f_post(4,i,1,k)
                fps_send_j(2,i,k) = f_post(9,i,1,k)
                fps_send_j(3,i,k) = f_post(10,i,1,k)
                fps_send_j(4,i,k) = f_post(16,i,1,k)
                fps_send_j(5,i,k) = f_post(18,i,1,k)
            enddo
        enddo

        call MPI_Bsend(fps_send_j, 5*nx*nz, MPI_DOUBLE_PRECISION, rank-nx_block, 3, MPI_COMM_WORLD, rc)
        call MPI_Recv(fps_recv_j, 5*nx*nz, MPI_DOUBLE_PRECISION, rank-nx_block, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
                
        ! scatter the received data (j=0)
        do k = 1, nz
            do i = 1, nx
                f_post(3,i,0,k) = fps_recv_j(1,i,k)
                f_post(7,i,0,k) = fps_recv_j(2,i,k)
                f_post(8,i,0,k) = fps_recv_j(3,i,k)
                f_post(15,i,0,k) = fps_recv_j(4,i,k)
                f_post(17,i,0,k) = fps_recv_j(5,i,k)
            enddo
        enddo
    endif

    ! /* exchange message with j=ny surface --- receive then send */
    if (block_y < ny_block) then
        ! collect data to send (j=ny)
        do k = 1, nz
            do i = 1, nx
                fps_send_j(1,i,k) = f_post(3,i,ny,k)
                fps_send_j(2,i,k) = f_post(7,i,ny,k)
                fps_send_j(3,i,k) = f_post(8,i,ny,k)
                fps_send_j(4,i,k) = f_post(15,i,ny,k)
                fps_send_j(5,i,k) = f_post(17,i,ny,k)
            enddo
        enddo

        call MPI_Bsend(fps_send_j, 5*nx*nz, MPI_DOUBLE_PRECISION, rank+nx_block, 4, MPI_COMM_WORLD, rc)
        call MPI_Recv(fps_recv_j, 5*nx*nz, MPI_DOUBLE_PRECISION, rank+nx_block, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
                
        ! scatter the received data (j=ny+1)
        do k = 1, nz
            do i = 1, nx
                f_post(4,i,ny+1,k) = fps_recv_j(1,i,k)
                f_post(9,i,ny+1,k) = fps_recv_j(2,i,k)
                f_post(10,i,ny+1,k) = fps_recv_j(3,i,k)
                f_post(16,i,ny+1,k) = fps_recv_j(4,i,k)
                f_post(18,i,ny+1,k) = fps_recv_j(5,i,k)
            enddo
        enddo
    endif

    ! /* exchange message with k=1 surface --- send then receive */
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

        call MPI_Bsend(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block, 5, MPI_COMM_WORLD, rc)
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

    ! /* exchange message with k=nz surface --- receive then send */
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

        call MPI_Bsend(fps_send_k, 5*nx*ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block, 6, MPI_COMM_WORLD, rc)
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


    !!! exchange messages with boundary lines -- message tag(for sender), velocity index in velocity model
    ! /* exchange with (nx, ny, k), message tag:7 */ 
    if (block_x < nx_block .AND. block_y < ny_block) then
        ! collect data to send (i=nx, j=ny)
        do k = 1, nz
            fpl_send_z(k) = f_post(7, nx, ny, k)
        enddo

        call MPI_Bsend(fpl_send_z, nz, MPI_DOUBLE_PRECISION, rank + nx_block + 1, 7, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_z, nz, MPI_DOUBLE_PRECISION, rank + nx_block + 1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)
        
        ! scatter the received data (i=nx+1, j=ny+1)
        do k = 1, nz
            f_post(10, nx+1, ny+1, k) = fpl_recv_z(k)
        enddo
    endif

    ! /* exchange with (1, ny, k), message tag:8 */ 
    if (block_x > 1 .AND. block_y < ny_block) then
        ! collect data to send (i=1, j=ny)
        do k = 1, nz
            fpl_send_z(k) = f_post(8, 1, ny, k)
        enddo

        call MPI_Bsend(fpl_send_z, nz, MPI_DOUBLE_PRECISION, rank + nx_block - 1, 8, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_z, nz, MPI_DOUBLE_PRECISION, rank + nx_block - 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=0, j=ny+1)
        do k = 1, nz
            f_post(9, 0, ny+1, k) = fpl_recv_z(k)
        enddo
    endif

    ! /* exchange with (nx, 1, k), message tag:9 */
    if (block_x < nx_block .AND. block_y > 1) then
        ! collect data to send (i=nx, j=1)
        do k = 1, nz
            fpl_send_z(k) = f_post(9, nx, 1, k)
        enddo

        call MPI_Bsend(fpl_send_z, nz, MPI_DOUBLE_PRECISION, rank - nx_block + 1, 9, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_z, nz, MPI_DOUBLE_PRECISION, rank - nx_block + 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=nx+1, j=0)
        do k = 1, nz
            f_post(8, nx+1, 0, k) = fpl_recv_z(k)
        enddo
    endif

    ! /* exchange with (1, 1, k), message tag:10 */
    if (block_x > 1 .AND. block_y > 1) then
        ! collect data to send (i=1, j=1)
        do k = 1, nz
            fpl_send_z(k) = f_post(10, 1, 1, k)
        enddo

        call MPI_Bsend(fpl_send_z, nz, MPI_DOUBLE_PRECISION, rank - nx_block - 1, 10, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_z, nz, MPI_DOUBLE_PRECISION, rank - nx_block - 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=0, j=0)
        do k = 1, nz
            f_post(7, 0, 0, k) = fpl_recv_z(k)
        enddo
    endif

    ! /* exchange with (nx, j, nz), message tag:11 */
    if (block_x < nx_block .AND. block_z < nz_block) then
        ! collect data to send (i=nx, k=nz)
        do j = 1, ny
            fpl_send_y(j) = f_post(11, nx, j, nz)
        enddo

        call MPI_Bsend(fpl_send_y, ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block + 1, 11, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_y, ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block + 1, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=nx+1, k=nz+1)
        do j = 1, ny
            f_post(14, nx+1, j, nz+1) = fpl_recv_y(j)
        enddo
    endif

    ! /* exchange with (1, j, nz), message tag:12 */
    if (block_x > 1 .AND. block_z < nz_block) then
        ! collect data to send (i=1, k=nz)
        do j = 1, ny
            fpl_send_y(j) = f_post(12, 1, j, nz)
        enddo

        call MPI_Bsend(fpl_send_y, ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block - 1, 12, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_y, ny, MPI_DOUBLE_PRECISION, rank + nx_block*ny_block - 1, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=0, k=nz+1)
        do j = 1, ny
            f_post(13, 0, j, nz+1) = fpl_recv_y(j)
        enddo
    endif

    ! /* exchange with (nx, j, 1), message tag:13 */
    if (block_x < nx_block .AND. block_z > 1) then
        ! collect data to send (i=nx, k=1)
        do j = 1, ny
            fpl_send_y(j) = f_post(13, nx, j, 1)
        enddo

        call MPI_Bsend(fpl_send_y, ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block + 1, 13, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_y, ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block + 1, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=nx+1, k=0)
        do j = 1, ny
            f_post(12, nx+1, j, 0) = fpl_recv_y(j)
        enddo
    endif

    ! /* exchange with (1, j, 1), message tag:14 */
    if (block_x > 1 .AND. block_z > 1) then
        ! collect data to send (i=1, k=1)
        do j = 1, ny
            fpl_send_y(j) = f_post(14, 1, j, 1)
        enddo

        call MPI_Bsend(fpl_send_y, ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block - 1, 14, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_y, ny, MPI_DOUBLE_PRECISION, rank - nx_block*ny_block - 1, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (i=0, k=0)
        do j = 1, ny
            f_post(11, 0, j, 0) = fpl_recv_y(j)
        enddo
    endif

    ! /* exchange with (i, ny, nz), message tag:15 */
    if (block_y < ny_block .AND. block_z < nz_block) then
        ! collect data to send (j=ny, k=nz)
        do i = 1, nx
            fpl_send_x(i) = f_post(15, i, ny, nz)
        enddo

        call MPI_Bsend(fpl_send_x, nx, MPI_DOUBLE_PRECISION, rank + nx_block + nx_block*ny_block, 15, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_x, nx, MPI_DOUBLE_PRECISION, rank + nx_block + nx_block*ny_block, 18, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (j=ny+1, k=nz+1)
        do i = 1, nx
            f_post(18, i, ny+1, nz+1) = fpl_recv_x(i)
        enddo
    endif
    
    ! /* exchange with (i, 1, nz), message tag:16 */
    if (block_y > 1 .AND. block_z < nz_block) then
        ! collect data to send (j=1, k=nz)
        do i = 1, nx
            fpl_send_x(i) = f_post(16, i, 1, nz)
        enddo

        call MPI_Bsend(fpl_send_x, nx, MPI_DOUBLE_PRECISION, rank - nx_block + nx_block*ny_block, 16, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_x, nx, MPI_DOUBLE_PRECISION, rank - nx_block + nx_block*ny_block, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (j=0, k=nz+1)
        do i = 1, nx
            f_post(17, i, 0, nz+1) = fpl_recv_x(i)
        enddo
    endif

    ! /* exchange with (i, ny, 1), message tag:17 */
    if (block_y < ny_block .AND. block_z > 1) then
        ! collect data to send (j=ny, k=1)
        do i = 1, nx
            fpl_send_x(i) = f_post(17, i, ny, 1)
        enddo

        call MPI_Bsend(fpl_send_x, nx, MPI_DOUBLE_PRECISION, rank + nx_block - nx_block*ny_block, 17, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_x, nx, MPI_DOUBLE_PRECISION, rank + nx_block - nx_block*ny_block, 16, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (j=ny+1, k=0)
        do i = 1, nx
            f_post(16, i, ny+1, 0) = fpl_recv_x(i)
        enddo
    endif

    ! /* exchange with (i, 1, 1), message tag:18 */
    if (block_y > 1 .AND. block_z > 1) then
        ! collect data to send (j=1, k=1)
        do i = 1, nx
            fpl_send_x(i) = f_post(18, i, 1, 1)
        enddo

        call MPI_Bsend(fpl_send_x, nx, MPI_DOUBLE_PRECISION, rank - nx_block - nx_block*ny_block, 18, MPI_COMM_WORLD, rc)
        call MPI_Recv(fpl_recv_x, nx, MPI_DOUBLE_PRECISION, rank - nx_block - nx_block*ny_block, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE, rc)

        ! scatter the received data (j=0, k=0)
        do i = 1, nx
            f_post(15, i, 0, 0) = fpl_recv_x(i)
        enddo
    endif

    call MPI_Buffer_detach(buffer, buffer_size, rc)
    
    deallocate(buffer)

end subroutine message_passing