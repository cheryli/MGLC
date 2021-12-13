!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

#define MRT
    module commondata
        implicit none
        
        integer, parameter :: nx=201,ny=nx
        integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
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


! #define BGK

    program main  
    use omp_lib   
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads

#ifdef BGK
    write(*,*) "--------------------------"
    write(*,*) "I am BGK collision operator"
    write(*,*) "--------------------------"
#endif
#ifdef MRT
    write(*,*) "--------------------------"
    write(*,*) "I am MRT collision operator"
    write(*,*) "--------------------------"
#endif

    call initial()
    !!call output_ASCII()

    !----------------------------------------OpenMP--------------------OpenMP----------
#ifdef _OPENMP
    call OMP_set_num_threads(2)
    myMaxThreads = OMP_get_max_threads()
    
    write(*,*) "Max Running threads=",myMaxThreads
    write(*,*) "    "
#endif
    !----------------------------------------OpenMP--------------------OpenMP----------
    
    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif
    
    do while((errorU.GT.eps).AND.(itc.LE.itc_max))

        itc = itc+1

        call collision()

        call streaming()

        call bounceback()

        call macro()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
        !if(MOD(itc,100000).EQ.0) then
        !    call output_Tecplot()
        !endif

    enddo
    
    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif
    
    write(*,*) "Time (CPU) = ", real(finish-start), "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", real(finish2-start2), "s"
#endif

    itc = itc+1
    call output_Tecplot()
    call output_binary()
    call getVelocity()
    
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

    !$omp parallel
    !$omp workshare
    rho = rho0
    u = 0.0d0
    v = 0.0d0
    !$omp end workshare
    !$omp end parallel
    
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
    
    !$omp parallel do default(none) shared(f,u,v,rho,ex,ey,omega) private(i,j,alpha,us2,un)
    do j=1,ny
        do i=1,nx
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine initial


    subroutine collision()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    !----------------------------
    real(8) :: us2, un(0:8)
    real(8) :: feq(0:8)
    !---------------------------
    real(8) :: s(0:8)
    real(8) :: m(0:8)
    real(8) :: m_post(0:8)
    real(8) :: meq(0:8)
    !---------------------------

#ifdef BGK
    !$omp parallel do default(none) shared(f,f_post,u,v,rho,ex,ey,omega) private(i,j,alpha,feq,us2,un)
    do j=1,ny
        do i=1,nx
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                feq(alpha) = omega(alpha)*rho(i,j)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                f_post(alpha,i,j) = f(alpha,i,j)-1.0d0/tauf*(f(alpha,i,j)-feq(alpha))
            enddo
        enddo
    enddo
    !$omp end parallel do
#endif

#ifdef MRT
    !$omp parallel do default(none) shared(f,f_post,rho,u,v) private(i,j,alpha,s,m,m_post,meq) 
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
    !$omp end parallel do
#endif

    return
    end subroutine collision


    subroutine streaming()
    use commondata
    implicit none
    integer :: i, j
    integer :: ip, jp
    integer :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey) private(i,j,ip,jp,alpha)
    do j=1,ny
        do i=1,nx
            do alpha=0,8
                ip = i-ex(alpha)
                jp = j-ey(alpha)
                
                f(alpha,i,j) = f_post(alpha,ip,jp)
                
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine streaming


    subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(f,f_post) private(j)
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
    !$omp end parallel do

    !$omp parallel do default(none) shared(f,f_post,Uwall,rho) private(i)
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
    !$omp end parallel do
    
    return
    end subroutine bounceback


    subroutine macro()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v) private(i,j)
    do j=1,ny
        do i=1,nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j) )/rho(i,j)
            v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j) )/rho(i,j)
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine macro


    subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2

    error1 = 0.0d0
    error2 = 0.0d0

    !$omp parallel do default(none) shared(u,up,v,vp) private(i,j) reduction(+:error1,error2)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
            
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j) 
        enddo
    enddo
    !$omp end parallel do

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

#ifdef BGK
    open(unit=02,file='BGKcavity-'//trim(filename)//'.dat',status='unknown')
#endif
#ifdef MRT
    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
#endif
    write(02,*) 'TITLE="Lid Driven Cavity"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" "Pressure" '
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xp(i), yp(j), u(i,j), v(i,j), rho(i,j)/3.0d0
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
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
#ifdef BGK
    open(41,file='BGKcavity-'//B2//'.plt',form='binary')
#endif
#ifdef MRT
    open(41,file='MRTcavity-'//B2//'.plt',form='binary')
#endif
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

