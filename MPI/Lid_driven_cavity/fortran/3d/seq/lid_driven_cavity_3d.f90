!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

module commondata
    implicit none
    
    integer, parameter :: nx=65, ny=65, nz =65
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
    ! call output_ASCII()
    call output_Tecplot()
    
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
            ! call getVelocity()
        endif

        if(MOD(itc,20000).EQ.0) then
            call getVelocity()
        endif

        if(MOD(itc,10000) .EQ. 0) then
            EXIT
        endif

    enddo
    
    call CPU_TIME(finish)
    
    write(*,*) "Time (CPU) = ", real(finish-start), "s"

    itc = itc+1
    call output_Tecplot()
    ! call output_ASCII()
    call output_binary()
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

    write(*,*) "nx=",nx,", ny=",ny, ",  nz = ",nz
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

    do j = 1, ny
        do i = 1, nx
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
    !---------------------------
    ! real(8) :: us2, un, feq
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
 
    f_post(0,i,j,k) = m_post(0)/19.0d0 - 5.0d0/399.0d0*m_post(1) + m_post(2)/21.0d0

    ! ------------
    f_post(1,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(3)/10.0d0 &
        - m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(3)/10.0d0 &
        + m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0
    
    ! ------------
    f_post(3,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(5)/10.0d0 &
        - m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(4,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(5)/10.0d0 &
        + m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(5,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(7)/10.0d0 &
        - m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(7)/10.0d0 &
        + m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    ! ------------
    f_post(7,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 + m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(8,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 - m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(9,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 + m_post(16)/8.0d0 + m_post(17)/8.0d0

    f_post(10,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 - m_post(16)/8.0d0 + m_post(17)/8.0d0

    ! -----------
    f_post(11,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 - m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(12,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 + m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(13,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 - m_post(16)/8.0d0 - m_post(18)/8.0d0
    
    f_post(14,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 + m_post(16)/8.0d0 - m_post(18)/8.0d0

    ! -----------
    f_post(15,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 + m_post(17)/8.0d0 - m_post(18)/8.0d0
    
    f_post(16,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 - m_post(17)/8.0d0 - m_post(18)/8.0d0

    f_post(17,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 + m_post(17)/8.0d0 + m_post(18)/8.0d0

    f_post(18,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 - m_post(17)/8.0d0 + m_post(18)/8.0d0


                ! us2 = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                ! do alpha = 0, 18
                !     un = u(i,j,k)*ex(alpha) + v(i,j,k)*ey(alpha) + w(i,j,k)*ez(alpha)
                !     feq = rho(i,j,k) * omega(alpha) * (1.0d0 + 3.0d0*un + 4.5d0*un*un - 1.5d0*us2)
                    
                !     f_post(alpha,i,j,k) = f(alpha,i,j,k) - Snu*(f(alpha,i,j,k) - feq)
                ! enddo
    
            enddo
        enddo
    enddo

    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp   
    integer :: alpha
    
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha=0,18
                    ip = i - ex(alpha)
                    jp = j - ey(alpha)
                    kp = k - ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

    do k = 1, nz
        do j = 1, ny
            ! (i = 1)
            f(1,1,j,k) = f_post(2,1,j,k)
            f(7,1,j,k) = f_post(10,1,j,k)
            f(9,1,j,k) = f_post(8,1,j,k)
            f(11,1,j,k) = f_post(14,1,j,k)
            f(13,1,j,k) = f_post(12,1,j,k)

            !(i = nx)
            f(2,nx,j,k) = f_post(1,nx,j,k)
            f(10,nx,j,k) = f_post(7,nx,j,k)
            f(8,nx,j,k) = f_post(9,nx,j,k)
            f(14,nx,j,k) = f_post(11,nx,j,k)
            f(12,nx,j,k) = f_post(13,nx,j,k)

        enddo
    enddo


    do k = 1, nz
        do i = 1, nx
            !(j = 1)
            f(3,i,1,k) = f_post(4,i,1,k)
            f(7,i,1,k) = f_post(10,i,1,k)
            f(8,i,1,k) = f_post(9,i,1,k)
            f(15,i,1,k) = f_post(18,i,1,k)
            f(17,i,1,k) = f_post(16,i,1,k)

            !(j = ny)
            f(4,i,ny,k) = f_post(3,i,ny,k)
            f(10,i,ny,k) = f_post(7,i,ny,k)
            f(9,i,ny,k) = f_post(8,i,ny,k)
            f(18,i,ny,k) = f_post(15,i,ny,k)
            f(16,i,ny,k) = f_post(17,i,ny,k)
        enddo
    enddo

    do j = 1, ny
        do i = 1, nx
            ! (k = 1)
            f(5,i,j,1) = f_post(6,i,j,1)
            f(11,i,j,1) = f_post(14,i,j,1)
            f(12,i,j,1) = f_post(13,i,j,1)
            f(15,i,j,1) = f_post(18,i,j,1)
            f(16,i,j,1) = f_post(17,i,j,1)

            ! (k = nz)
            f(6,i,j,nz) = f_post(5,i,j,nz)
            f(14,i,j,nz) = f_post(11,i,j,nz) - rho(i,j,nz) / 6.0d0 * (U0)
            f(13,i,j,nz) = f_post(12,i,j,nz) - rho(i,j,nz) / 6.0d0 * (-U0)
            f(18,i,j,nz) = f_post(15,i,j,nz)
            f(17,i,j,nz) = f_post(16,i,j,nz)
        enddo
    enddo
    
    return
end subroutine bounceback


subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k, alpha

    rho = 0.0d0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 18
                    rho(i,j,k) = rho(i,j,k) + f(alpha,i,j,k)
                    u(i,j,k) = u(i,j,k) + f(alpha,i,j,k) * ex(alpha)
                    v(i,j,k) = v(i,j,k) + f(alpha,i,j,k) * ey(alpha)
                    w(i,j,k) = w(i,j,k) + f(alpha,i,j,k) * ez(alpha)
                enddo

                u(i,j,k) = u(i,j,k)/rho(i,j,k)
                v(i,j,k) = v(i,j,k)/rho(i,j,k)
                w(i,j,k) = w(i,j,k)/rho(i,j,k)
            enddo
        enddo
    enddo

    return
end subroutine macro


subroutine check()
    use commondata
    implicit none
    integer :: i, j, k
    real(8) :: error1, error2

    error1 = 0.0d0
    error2 = 0.0d0

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                error1 = error1 + (u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k))+(v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k))
                error2 = error2 + u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
            enddo
        enddo
    enddo

    up = u
    vp = v
    wp = w

    errorU = sqrt(error1)/sqrt(error2)

    write(*,*) itc,' ',errorU

    return
end subroutine check


subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j, k
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=77,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
        write(77,*) 'TITLE="Lid Driven Cavity"'
        write(77,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "Pressure" '
        write(77,101) nx, ny, nz
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(77,100) xp(i), yp(j), zp(k), u(i,j,k), v(i,j,k), w(i,j,k), rho(i,j,k)/3.0d0
                enddo
            enddo
        enddo
100     format(1x,3(e11.4,' '),10(e13.6,' '))
101     format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'K=',1x,i5,1x,'F=POINT')
    close(77)

    return
end subroutine output_ASCII


subroutine output_binary()
    use commondata
    implicit none
    integer :: i, j, k
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)
    
    open(unit=01,file='MRTcavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((v(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((rho(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
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
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    
    write(B2,'(i9.9)') itc

    open(41,file='MRTcavity-'//B2//'.plt', access='stream', form='unformatted')
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
    write(41) 7

    !c-- Variable names.
    V1 = 'X'
    call dumpstring(V1)
    V2 = 'Y'
    call dumpstring(V2)
    V3 = 'Z'
    call dumpstring(V3)
    V4 = 'U'
    call dumpstring(V4)
    V5 = 'V'
    call dumpstring(V5)
    V6 = 'W'
    call dumpstring(V6)
    V7 = 'Pressure'
    call dumpstring(V7)

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
    write(41) nz

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

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(zp(k))
                write(41) real(u(i,j,k))
                write(41) real(v(i,j,k))
                write(41) real(w(i,j,k))
                write(41) real(rho(i,j,k)/3.0d0)
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
    integer :: i, j, k
    character(len=9) :: B2

    write(B2,'(i9.9)') itc

    open(unit=02,file='./u-z_'//B2//'.dat',status='unknown')
    do k=1,nz
        write(02,*) u(nxHalf, nyHalf, k)/U0, zp(k)/dble(nz)
    enddo
    close(02)

    open(unit=03,file='./x-w_'//B2//'.dat',status='unknown')
    do i=1,nx
        write(03,*) xp(i)/dble(nx), w(i, nyHalf, nzHalf)/U0
    enddo
    close(03)

    return
end subroutine getVelocity

