module commondata
    implicit none
    
    integer, parameter :: total_nx=121, total_ny=121, total_nz = 121
    integer, parameter :: nx_block = 2, ny_block = 2, nz_block = 2
    integer :: nx, ny, nz  
    integer :: block_x, block_y, block_z
    real(8), parameter :: Reynolds=1000.0d0
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: U0=0.1d0
    real(8), parameter :: tauf=U0*dble(total_nx)/Reynolds*3.0d0+0.5d0
    
    integer :: itc
    integer, parameter :: itc_max=INT(50000000)
    
    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    
    real(8) :: xp(0:total_nx+1), yp(0:total_ny+1), zp(0:total_nz+1)
    real(8), allocatable :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(8), allocatable :: rho(:, :, :)
    real(8), allocatable :: up(:, :, :), vp(:, :, :), wp(:, :, :)
    
    real(8), allocatable :: f(:, :, :, :), f_post(:, :, :, :)
    
    integer :: idx
    ! real(8), parameter :: omega(0:18) = (/ 1.0d0/3.0d0, & 
    !         1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, &
    !         1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
    !         1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
    real(8), parameter :: omega(0:18) = (/ 1.0d0/3.0d0, &
                                        (1.0d0/18.0d0, idx = 1, 6), &
                                        (1.0d0/36.0d0, idx = 7, 18) /) 
    integer, parameter :: ex(0:18) = (/ 0,  &
                                        1, -1,  0,  0,  0,  0, &
                                        1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
    integer, parameter :: ey(0:18) = (/ 0, &
                                        0,  0,  1, -1,  0,  0, &
                                        1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)  
    integer, parameter :: ez(0:18) = (/ 0, &
                                        0,  0,  0,  0,  1, -1, &
                                        0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    integer :: rc, rank, num_process


    
end module commondata