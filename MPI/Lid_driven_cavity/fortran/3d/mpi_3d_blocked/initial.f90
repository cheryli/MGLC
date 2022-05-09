subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(8) :: un(0:18)
    real(8) :: us2

    itc = 0
    errorU = 100.0d0

    if (rank == 0) then
        write(*,*) "nx=",nx,", ny=",ny, ",  nz = ",nz
        write(*,*) "Reynolds=",real(Reynolds)
        write(*,*) "U0=",real(U0),",    tauf=",real(tauf)
    endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i=1,total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j=1,total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(total_nz+1) = dble(total_nz)
    do k=1,total_nz
        zp(k) = dble(k)-0.5d0
    enddo

    
    allocate (u(nx, ny, nz))
    allocate (v(nx, ny, nz))
    allocate (w(nx, ny, nz))
    allocate (rho(nx, ny, nz))
    allocate (up(nx, ny, nz))
    allocate (vp(nx, ny, nz))
    allocate (wp(nx, ny, nz))
    
    allocate (f(0:18, nx, ny, nz))
    allocate (f_post(0:18, 0:nx+1, 0:ny+1, 0:nz+1))

    rho = rho0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    up = 0.0d0
    vp = 0.0d0
    wp = 0.0d0

    ! top bpundary
    if (coords(2) == dims(2) - 1) then
        do j = 1, ny
            do i = 1, nx
                u(i, j, nz) = U0
            enddo
        enddo
    endif

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                us2 = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                do alpha = 0, 18
                    un(alpha) = u(i,j,k)*ex(alpha) + v(i,j,k)*ey(alpha) + w(i,j,k)*ez(alpha)
                    f(alpha,i,j,k) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine initial