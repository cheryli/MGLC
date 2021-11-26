!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

module commondata
    implicit none
    
    integer, parameter :: nx=201, ny=201, nz = 201
    integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1, nzHalf=(nz-1)/2+1
    real(8), parameter :: Reynolds=1000.0d0
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: U0=0.1d0
    real(8), parameter :: tauf=U0*dble(nx)/Reynolds*3.0d0+0.5d0
    
    integer :: itc
    integer, parameter :: itc_max=INT(50000000)
    
    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    
    real(8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
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
end module commondata



program main
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads

    write(*,*) "--------------------------"
    write(*,*) "I am MRT collision operator"
    write(*,*) "--------------------------"


    call initial()
    
    call CPU_TIME(start)
    
    do while((errorU.GT.eps).AND.(itc.LE.itc_max))

        itc = itc+1

        call collision()

        call streaming()

        call bounceback()

        call macro()

        if(MOD(itc,2000).EQ.0) then
            call check()
            ! call output_ASCII()
        endif
        !if(MOD(itc,100000).EQ.0) then
        !    call output_Tecplot()
        !endif

    enddo
    
    call CPU_TIME(finish)
    
    write(*,*) "Time (CPU) = ", real(finish-start), "s"

    itc = itc+1
    call output_Tecplot()
    call output_ASCII()
    ! call output_binary()
    call getVelocity()
    
    write(*,*) "Deallocate Array..."
    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(wp)
    deallocate(f)
    deallocate(f_post)
    write(*,*) "    "
    
    write(*,*) "Successfully: DNS completed!"

    stop
end program main


subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(8) :: un(0:18)
    real(8) :: us2

    itc = 0
    errorU = 100.0d0

    write(*,*) "nx=",nx,", ny=",ny
    write(*,*) "Reynolds=",real(Reynolds)
    write(*,*) "U0=",real(U0),",    tauf=",real(tauf)

    xp(0) = 0.0d0
    xp(nx+1) = dble(nx)
    do i=1,nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(ny+1) = dble(ny)
    do j=1,ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(nz+1) = dble(nz)
    do k=1,nz
        yp(k) = dble(k)-0.5d0
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

    do j = 1, ny
        do i=1, nx
            u(i, j, nz) = U0
        enddo
    enddo

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


subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !---------------------------
    real(8) :: s(0:18)
    real(8) :: m(0:18)
    real(8) :: m_post(0:18)
    real(8) :: meq(0:18)
    real(8) :: us2
    !---------------------------
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
    ! loop body, re-indented for readability 
        
    m(0) = f(0,i,j,k) &
        + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(1) = -30.0d0 * f(0,i,j,k) &
        - 11.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + 8.0d0 * (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k)) &
        + 8.0d0 * (f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(2) = 12.0d0 * f(0,i,j,k) &
        - 4.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(3) = f(1,i,j,k) - f(2,i,j,k) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)
    m(4) = -4.0d0 * (f(1,i,j,k) - f(2,i,j,k)) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)       
    m(5) = f(3,i,j,k) - f(4,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(6) = -4.0d0 * (f(3,i,j,k) - f(4,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(7) = f(5,i,j,k) - f(6,i,j,k) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(8) = -4.0d0 * (f(5,i,j,k) - f(6,i,j,k)) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(9) = 2.0d0 * (f(1,i,j,k) + f(2,i,j,k)) - f(3,i,j,k) - f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(10) = -4.0d0 * (f(1,i,j,k) + f(2,i,j,k)) + 2.0d0 * (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(11) = f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(12) = -2.0d0 * (f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(13) = f(7,i,j,k) - f(8,i,j,k) - f(9,i,j,k) + f(10,i,j,k)
    m(14) = f(15,i,j,k) - f(16,i,j,k) - f(17,i,j,k) + f(18,i,j,k)
    m(15) = f(11,i,j,k) - f(12,i,j,k) - f(13,i,j,k) + f(14,i,j,k)
    m(16) = f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) - f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) + f(14,i,j,k)
    m(17) = -f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(18) = f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) - f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
        
        
        meq(0) = rho(i,j,k) 
        meq(1) = rho(i,j,k) * (-11.0d0 + 19.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(2) = rho(i,j,k) * (3.0d0 - 11.0d0/2.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(3) = rho(i,j,k) * u(i,j,k)
        meq(4) = -2.0d0/3.0d0 * rho(i,j,k) * u(i,j,k)
        meq(5) = rho(i,j,k) * v(i,j,k)
        meq(6) = -2.0d0/3.0d0 * rho(i,j,k) * v(i,j,k)
        meq(7) = rho(i,j,k) * w(i,j,k)
        meq(8) = -2.0d0/3.0d0 * rho(i,j,k) * w(i,j,k)
        meq(9) = rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(10) = -1.0d0/2.0d0 * rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(11) = rho(i,j,k) * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(12) = -1.0d0/2.0d0 * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(13) = rho(i,j,k) * u(i,j,k) * v(i,j,k)
        meq(14) = rho(i,j,k) * v(i,j,k) * w(i,j,k)
        meq(15) = rho(i,j,k) * u(i,j,k) * w(i,j,k)
        meq(16) = 0.0d0
        meq(17) = 0.0d0
        meq(18) = 0.0d0


        s(0) = 0.0d0    !! s_{\rho}
        s(1) = Snu      !! s_{e}
        s(2) = Snu      !! s_{epsilon}
        s(3) = 0.0d0    !! s_{j}
        s(4) = Sq       !! s_{q}
        s(5) = 0.0d0    !! s_{j}
        s(6) = Sq       !! s_{q}
        s(7) = 0.0d0    !! s_{j}
        s(8) = Sq       !! s_{q}
        s(9) = Snu      !! s_{\nu}
        s(10) = Snu     !! s_{\pi}
        s(11) = Snu     !! s_{\nu}
        s(12) = Snu     !! s_{\pi}
        s(13) = Snu     !! s_{nu}
        s(14) = Snu     !! s_{nu}
        s(15) = Snu     !! s_{nu}
        s(16) = Sq      !! s_{m}
        s(17) = Sq      !! s_{m}
        s(18) = Sq      !! s_{m}
        
        do alpha=0,18
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
        enddo

    f_post(0,i,j,k) = m_post(0)/19.0d0 - 5.0d0/399.0d0*
    f_post(1,i,j,k) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j,k) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0 
    f_post(3,i,j,k) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j,k) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j,k) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j,k) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j,k) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j,k) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
            enddo
        enddo
    enddo

    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j
    integer :: ip, jp
    integer :: alpha
    
    do j=1,ny
        do i=1,nx
            do alpha=0,8
                ip = i-ex(alpha)
                jp = j-ey(alpha)
                
                f(alpha,i,j) = f_post(alpha,ip,jp)
                
            enddo
        enddo
    enddo
    
    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j

    do j=1,ny
        !Left side (i=1)
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side (i=nx)
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo


    do i=1,nx
        !Bottom side (j=1)
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)
        f(6,i,1) = f_post(8,i,1)

        !Top side (j=ny)
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*(U0)/6.0d0
        f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-U0)/6.0d0
    enddo
    
    return
end subroutine bounceback


subroutine macro()
    use commondata
    implicit none
    integer :: i, j

    do j=1,ny
        do i=1,nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j) )/rho(i,j)
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j) )/rho(i,j)
        enddo
    enddo

    return
end subroutine macro


subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2

    error1 = 0.0d0
    error2 = 0.0d0

    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
            
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j) 
        enddo
    enddo

    errorU = sqrt(error1)/sqrt(error2)

    write(*,*) itc,' ',errorU

    return
end subroutine check


subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
        write(02,*) 'TITLE="Lid Driven Cavity"'
        write(02,*) 'VARIABLES="X" "Y" "U" "V" "Pressure" '
        write(02,101) nx, ny
        do j=1,ny
            do i=1,nx
                write(02,100) xp(i), yp(j), u(i,j), v(i,j), rho(i,j)/3.0d0
            enddo
        enddo
100     format(1x,2(e11.4,' '),10(e13.6,' '))
101     format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII


subroutine output_binary()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)
    
    open(unit=01,file='MRTcavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) ((u(i,j),i=1,nx),j=1,ny)
    write(01) ((v(i,j),i=1,nx),j=1,ny)
    write(01) ((rho(i,j),i=1,nx),j=1,ny)
    close(01)

    return
    end subroutine output_binary


!!!c--------------------------------
!!!c--------------------------------
    subroutine output_Tecplot()
    use commondata
    implicit none
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5
    integer, parameter :: kmax=1
    character(len=40) :: zoneName
    
    write(B2,'(i9.9)') itc

    open(41,file='MRTcavity-'//B2//'.plt',form='binary')
    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"

    !c--Integer value of 1
    write(41) 1

    Title="MyFirst"
    call dumpstring(title)

    !c-- Number of variables in this data file (here 5 variables)
    write(41) 5

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='Pressure'
    call dumpstring(V5)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName="ZONE 001"
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
                write(41) real(rho(i,j)/3.0d0)
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
    end subroutine output_Tecplot
!!!c--------------------------------
    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer :: stringLength
    integer :: ii
    integer :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring
!!!c--------------------------------
!!!c-------------------------------- 


subroutine getVelocity()
    use commondata
    implicit none
    integer :: i, j

    open(unit=02,file='./u-y.dat',status='unknown')
    do j=1,ny
        write(02,*) u(nxHalf,j)/U0, yp(j)/dble(nx)
    enddo
    close(02)
    
    open(unit=03,file='./x-v.dat',status='unknown')
    do i=1,nx
        write(03,*) xp(i)/dble(nx), v(i,nyHalf)/U0
    enddo
    close(03)

    return
end subroutine getVelocity

