subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

    if (block_x == 1) then
        do k = 1, nz
            do j = 1, ny
                ! (i = 1)
                f(1,1,j,k) = f_post(2,1,j,k)
                f(7,1,j,k) = f_post(10,1,j,k)
                f(9,1,j,k) = f_post(8,1,j,k)
                f(11,1,j,k) = f_post(14,1,j,k)
                f(13,1,j,k) = f_post(12,1,j,k)
            enddo
        enddo
    endif

    if (block_x == nx_block) then
        do k = 1, nz
            do j = 1, ny
                !(i = nx)
                f(2,nx,j,k) = f_post(1,nx,j,k)
                f(10,nx,j,k) = f_post(7,nx,j,k)
                f(8,nx,j,k) = f_post(9,nx,j,k)
                f(14,nx,j,k) = f_post(11,nx,j,k)
                f(12,nx,j,k) = f_post(13,nx,j,k)
            enddo
        enddo
    endif


    if (block_y == 1) then
        do k = 1, nz
            do i = 1, nx
                !(j = 1)
                f(3,i,1,k) = f_post(4,i,1,k)
                f(7,i,1,k) = f_post(10,i,1,k)
                f(8,i,1,k) = f_post(9,i,1,k)
                f(15,i,1,k) = f_post(18,i,1,k)
                f(17,i,1,k) = f_post(16,i,1,k)
            enddo
        enddo
    endif

    if (block_y == ny_block) then
        do k = 1, nz
            do i = 1, nx
                !(j = ny)
                f(4,i,ny,k) = f_post(3,i,ny,k)
                f(10,i,ny,k) = f_post(7,i,ny,k)
                f(9,i,ny,k) = f_post(8,i,ny,k)
                f(18,i,ny,k) = f_post(15,i,ny,k)
                f(16,i,ny,k) = f_post(17,i,ny,k)
            enddo
        enddo
    endif

    if (block_z == 1) then
        do j = 1, ny
            do i = 1, nx
                ! (k = 1)
                f(5,i,j,1) = f_post(6,i,j,1)
                f(11,i,j,1) = f_post(14,i,j,1)
                f(12,i,j,1) = f_post(13,i,j,1)
                f(15,i,j,1) = f_post(18,i,j,1)
                f(16,i,j,1) = f_post(17,i,j,1)
            enddo
        enddo
    endif

    if (block_z == nz_block) then
        do j = 1, ny
            do i = 1, nx
                ! (k = nz)
                f(6,i,j,nz) = f_post(5,i,j,nz)
                f(14,i,j,nz) = f_post(11,i,j,nz) - rho(i,j,nz) / 6.0d0 * (U0)
                f(13,i,j,nz) = f_post(12,i,j,nz) - rho(i,j,nz) / 6.0d0 * (-U0)
                f(18,i,j,nz) = f_post(15,i,j,nz)
                f(17,i,j,nz) = f_post(16,i,j,nz)
            enddo
        enddo
    endif
    
    return
end subroutine bounceback