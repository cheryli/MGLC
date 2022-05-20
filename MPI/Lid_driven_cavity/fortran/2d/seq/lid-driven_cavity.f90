!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

module commondata
    implicit none
    
    integer, parameter :: nx=201,ny=201
    real(8), parameter :: Reynolds=1000.0d0
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: U0=0.1d0
    real(8), parameter :: tauf=U0*dble(nx)/Reynolds*3.0d0+0.5d0
    
    integer :: itc
    integer, parameter :: itc_max=INT(50000000)
    
    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    
    real(8) :: xp(0:nx+1), yp(0:ny+1)
    real(8), allocatable :: u(:,:), v(:,:)
    real(8), allocatable :: rho(:,:)
    real(8), allocatable :: up(:,:), vp(:,:)
    
    real(8), allocatable :: f(:,:,:), f_post(:,:,:)
    
    real(8) :: omega(0:8)
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
    real(8) :: Uwall(1:2)
    data Uwall/U0,  0/
    
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

    call output()
    
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

    call output()
    
    write(*,*) "Deallocate Array..."
    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(f)
    deallocate(f_post)
    write(*,*) "    "
    
    write(*,*) "Successfully: DNS completed!"

    stop
end program main


subroutine initial()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
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
    
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (rho(nx,ny))
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    
    allocate (f(0:8,nx,ny))
    allocate (f_post(0:8,0:nx+1,0:ny+1))

    rho = rho0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0
    
    do i=1,nx
        u(i,ny) = U0
    enddo

    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    do j=1,ny
        do i=1,nx
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
            enddo
        enddo
    enddo

    return
end subroutine initial


subroutine collision()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    !---------------------------
    real(8) :: s(0:8)
    real(8) :: m(0:8)
    real(8) :: m_post(0:8)
    real(8) :: meq(0:8)
    !---------------------------

    do j=1,ny
        do i=1,nx
        
    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
    m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
            meq(8) = rho(i,j)*( u(i,j)*v(i,j) )

            s(0) = 0.0d0      !!s_{\rho}
            s(1) = Snu !!s_{e}
            s(2) = Snu !!s_{\epsilon}
            s(3) = 0.0d0      !!s_{j} 
            s(4) = Sq !!s_{q}
            s(5) = 0.0d0      !!s_{j}
            s(6) = Sq       !!s_{q}
            s(7) = Snu !!s_{\nu}
            s(8) = Snu       !!s_{\nu}
            
            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0 
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
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
        f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*(Uwall(1))/6.0d0
        f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-Uwall(1))/6.0d0
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

subroutine output()
    use commondata
    implicit none
    real(8) :: stream(nx, ny), vorticity(nx, ny)

    call compute_stream_vorticity(stream, vorticity, u, v, nx, ny, U0)

    ! call output_ASCII(xp, yp, u, v, rho, nx, ny, itc)
    call output_binary(u, v, rho, nx, ny, itc)
    call output_Tecplot(xp, yp, u, v, rho, stream, vorticity, nx, ny, itc)
    call getVelocity(xp, yp, u, v, rho, nx, ny, U0, itc)

    return
end subroutine output

subroutine compute_stream_vorticity(stream, vorticity, u, v, nx, ny, U0)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: u(nx, ny), v(nx, ny), U0
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
        vorticity(i, ny) = 2.0d0 * (U0 - u(i, ny-1))
    enddo

    vorticity(1, ny) = 2.0d0 * U0
    vorticity(nx, ny) = 2.0d0 * U0

    vorticity(1, 1) = 0
    vorticity(nx, 1) = 0

end subroutine compute_stream_vorticity



subroutine output_ASCII(xp, yp, u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
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


subroutine output_binary(u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
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
subroutine output_Tecplot(xp, yp, u, v, rho, stream, vorticity, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), stream(nx, ny), vorticity(nx, ny)
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1, V2, V3, V4, V5, V6, V7
    integer, parameter :: kmax=1
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
    V6 = 'Stream'
    call dumpstring(V6)
    V7 = 'Vorticity'
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
                write(41) real(stream(i, j)) 
                write(41) real(vorticity(i, j))
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


subroutine getVelocity(xp, yp, u, v, rho, nx, ny, U0, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), U0
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer :: i, j
    integer :: nxHalf, nyHalf

    nxHalf = (nx-1)/2 + 1
    nyHalf = (ny-1)/2 + 1

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

