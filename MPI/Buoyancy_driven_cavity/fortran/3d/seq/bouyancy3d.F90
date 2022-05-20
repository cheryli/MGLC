!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

    module commondata
        implicit none
        
        integer, parameter :: loadInitField=0
        integer, parameter :: flowReversal=0  !! value 0: do NOT simulate reversal; value 1: do simulate reversal
    
        integer, parameter :: nx=51, ny=nx, nz=nx
        integer, parameter :: nxHalf=(nx-1)/2+1,nyHalf=(ny-1)/2+1,nzHalf=(nz-1)/2+1
        
        real(kind=8), parameter :: Rayleigh=1e6
        real(kind=8), parameter :: Prandtl=0.71d0
        real(kind=8), parameter :: Mach=0.1d0

        real(kind=8), parameter :: tauf=0.5d0+Mach*dble(nz)*DSQRT(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        
        real(kind=8), parameter :: Ekman=0.001d0
        real(kind=8), parameter :: omegaRatating=viscosity/2.0d0/Ekman/dble(nz*nz)
        
        real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/dble(nz)
        real(kind=8), parameter :: gBeta=gBeta1/dble(nz*nz)
        
        real(kind=8), parameter :: timeUnit=dsqrt(dble(nz)/gBeta)  !!dble(ny*ny)/diffusivity
        integer, parameter :: dimensionlessTimeMax=1000
        integer, parameter :: flowReversalTime=20000
        integer :: itc
        integer, parameter :: itc_max=INT(dimensionlessTimeMax*timeUnit)
        
        real(kind=8), parameter :: epsU=1e-8
        real(kind=8), parameter :: epsT=1e-8
        real(kind=8) :: errorU, errorT
        
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:)
        real(kind=8), allocatable :: rho(:,:,:)
        real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)

        real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
        real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
        
        integer :: ex(0:18), ey(0:18), ez(0:18)
        data ex/0, 1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0/
        data ey/0, 0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1/
        data ez/0, 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1/
        
        real(kind=8), allocatable :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
        
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
        
        integer :: fileNum
        integer :: dimensionlessTime
        integer :: statisticallyStationaryState, statisticallyStationaryTime
        integer :: statisticallyConvergeState
        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
        real(kind=8) :: NuVolAvg_mean(0:dimensionlessTimeMax), ReVolAvg_mean(0:dimensionlessTimeMax)
        
    end module commondata
        
!-----------------------------------------------------------------------------
#define outputTecplot
#define outputBinFile

!!==velocity B.C.==
#define noslipWalls
!!==velocity B.C.==

!~~temperature B.C. (for RB convection)~~
!#define RBconvection
!#define BackFrontWallsAdiabatic
!#define LeftRightWallsAdiabatic
!#define TopBottomPlatesConstT
!~~temperature B.C.~~

!~~temperature B.C. (for cavity flow benchmark)~~
#define benchmarkCavity
#define BackFrontWallsAdiabatic
#define LeftRightWallsConstT
#define TopBottomPlatesAdiabatic
!~~temperature B.C.~~

    program main
    use omp_lib    
    use commondata
    implicit none
    real(kind=8) :: start, finish
    real(kind=8) :: start2, finish2
    integer :: myMaxThreads

#ifdef _OPENMP
    write(*,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(24)
    myMaxThreads = OMP_get_max_threads()
    write(*,100) "|----------Max Running threads:",myMaxThreads,"----------|"
100 format(1X,A,I5)
#endif


    call initial()
#ifdef outputTecplot
    call output_Tecplot()
#endif

    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max).AND.(statisticallyConvergeState.EQ.0) )

        itc = itc+1

        call collision()

        call streaming()

        call bounceback()

        call collisionT()

        call streamingT()

        call bouncebackT()

        call macro()

        call macroT()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

#ifdef benchmarkCavity
        if(MOD(itc,1000000).EQ.0) then
#ifdef outputBinFile
            call backupData()
#endif
        endif
#endif
        
#ifdef RBconvection
        if( MOD(itc,int(timeUnit)).EQ.0 ) then
            
            call calNuRe()
            
            if(statisticallyStationaryState.EQ.1) then
#ifdef outputBinFile
                call output_binary()
#endif
            endif
            
        endif
#endif

    enddo
    
    write(*,*) "Backup Distribution Function......"
#ifdef outputBinFile
    call backupData()
#endif
#ifdef outputTecplot
            call output_Tecplot()
#endif
    write(*,*) "    "
    
    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif
    write(*,*) "---------------------------------------------"
    write(*,*) "Time (CPU) = ", real(finish-start), "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", real(finish2-start2), "s"
#endif
    write(*,*) "---------------------------------------------"

    write(*,*) "Dellocate Array......"
    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(T)
    deallocate(up)
    deallocate(vp)
    deallocate(wp)
    deallocate(Tp)
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Fz)
    write(*,*) "    "
    
    write(*,*) "Successfully: DNS completed!"
    end program main
    
    
    subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: un(0:18), unT(0:6)
    real(kind=8) :: us2
    real(kind=8) :: omega(0:18), omegaT(0:6)
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    
        write(*,*) 'Mesh:',nx,ny,nz
        write(*,*) 'Rayleigh=',real(Rayleigh), ', Prandtl=',real(Prandtl), ', Mach=',real(Mach)
        write(*,*) "Ekman=",real(Ekman)
        write(*,*) "   "
        write(*,*) 'tauf=',real(tauf)
        write(*,*) "paraA=",real(paraA)
        write(*,*) "viscosity=",real(viscosity), ", diffusivity=",real(diffusivity)
        write(*,*) "omegaRatating=", real(omegaRatating)
        write(*,*) "itc_max=",itc_max
        write(*,*) "Ouput will begin at", int(timeUnit)
        write(*,*) "Output interval is", int(timeUnit)
        write(*,*) "Time unit: Sqrt(L0/(gBeta*DeltaT))=",dsqrt(dble(Ny)/gBeta)
        write(*,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT)=",dsqrt(gBeta*dble(Ny))
        write(*,*) "    "
    
#ifdef RBconvection
    write(*,*) "    "
    write(*,*) "I am RBconvection"
    if(flowReversal.EQ.1) then
        write(*,*) "I am also flow reversal"
    write(*,*) "    "
    elseif(flowReversal.EQ.0) then
        write(*,*) "I am NOT flow reversal"
    write(*,*) "    "
    else
        write(*,*) "Error: Please check flowReversal setting!"
        stop
    endif
#endif

#ifdef benchmarkCavity
    write(*,*) "I am benchmarkCavity"
    write(*,*) "    "
#endif

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
    
    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (T(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    allocate (up(nx,ny,nz))
    allocate (vp(nx,ny,nz))
    allocate (wp(nx,ny,nz))
    allocate (Tp(nx,ny,nz))
    
    allocate (f(0:18,nx,ny,nz))
    allocate (f_post(0:18,0:nx+1,0:ny+1,0:nz+1))
    allocate (g(0:6,nx,ny,nz))
    allocate (g_post(0:6,0:nx+1,0:ny+1,0:nz+1))
    
    allocate (Fx(nx,ny,nz))
    allocate (Fy(nx,ny,nz))
    allocate (Fz(nx,ny,nz))
    
    rho = 1.0d0
    
    omega(0) = 1.0d0/3.0d0
    do alpha=1,6
        omega(alpha) = 1.0d0/18.0d0
    enddo
    do alpha=7,18
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/7.0d0
    do alpha=1,6
        omegaT(alpha) = (paraA+6.0d0)/42.0d0
    enddo

    if(loadInitField.EQ.0) then 
    
        write(*,*) "Initial field is set exactly"
        u = 0.0d0
        v = 0.0d0
        w = 0.0d0
        T = 0.0d0
        
#ifdef noslipWalls
    write(*,*) "Velocity B.C. for left/right walls are: ===No-slip wall==="
#endif
    
#ifdef LeftRightWallsConstT
    do k=1,nz
        do i=1,nx
            T(i,1,k) = Thot
            T(i,ny,k) = Tcold
        enddo
    enddo
    write(*,*) "Temperature B.C. for left wall is:===Hot left wall==="
    write(*,*) "Temperature B.C. for right wall is:==Cold right wall==="
#endif
#ifdef LeftRightWallsAdiabatic
    write(*,*) "Temperature B.C. for left/right walls are:===Adiabatic walls==="
#endif


#ifdef TopBottomPlatesConstT
    do j=1,ny
        do i=1,nx
            T(i,j,1) = Thot
            T(i,j,nz) = Tcold
        enddo
    enddo
    write(*,*) "Temperature B.C. for bottom plate is:==Hot bottom plate==="
    write(*,*) "Temperature B.C. for top plate is:====Cold top plate==="
#endif
#ifdef TopBottomPlatesAdiabatic
    write(*,*) "Temperature B.C. for top/bottom plates are:===Adiabatic plates==="
#endif

#ifdef BackFrontWallsAdiabatic
    write(*,*) "Temperature B.C. for back/front walls are:===Adiabatic walls==="
#endif
        
        !$omp parallel do default(none) shared(f,g,u,v,w,T,ex,ey,ez,omega,omegaT,rho) private(i,j,k,alpha,us2,un,unT)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                    do alpha=0,18
                        un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        f(alpha,i,j,k) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                    enddo
                    do alpha=0,6
                        unT(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        g(alpha,i,j,k) = omegaT(alpha)*T(i,j,k)*(1.0d0+21.0d0/(6.0d0+paraA)*unT(alpha))
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

    elseif(loadInitField.EQ.1) then
        write(*,*) "Load initial field from previous simulation"
        write(*,*) "Start to read raw data >>>>>>>>>>>>"
        open(unit=01,file='./reload/backupFile-1000.bin',form="unformatted",access="sequential",status='old')
        read(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(01) ((((f(alpha,i,j,k),alpha=0,18),i=1,nx),j=1,ny),k=1,nz)
        read(01) ((((g(alpha,i,j,k),alpha=0,6),i=1,nx),j=1,ny),k=1,nz)
        close(01)
        write(*,*) "Read raw data!"
    else
        write(*,*) "Error: initial field is not properly set"
    endif
    
    up = 0.0d0
    vp = 0.0d0
    wp = 0.0d0
    Tp = 0.0d0
    
    f_post = 0.0d0
    g_post = 0.0d0
        
    fileNum = 0
    dimensionlessTime = 0
    
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 0
    statisticallyConvergeState = 0

    return
    end subroutine initial
    
    
    subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,w,Fx,Fy,Fz,T) private(i,j,k,alpha,s,m,m_post,meq,fSource) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
    !--------------------------------------------------------------------------------------------------------------------
    !---m0    
    m(0) =f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m1
    m(1) = -30.0d0*f(0,i,j,k)-11.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
    +8.0d0*( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m2
    m(2) = 12.0d0*f(0,i,j,k)-4.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m3
    m(3) = f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m4
    m(4) = -4.0d0*(f(1,i,j,k)-f(2,i,j,k))+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k) &
                    +f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m5
    m(5) = f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m6
    m(6) = -4.0d0*(f(3,i,j,k)-f(4,i,j,k))+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k) &
                +f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m7
    m(7) = f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m8
    m(8) = -4.0d0*(f(5,i,j,k)-f(6,i,j,k))+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k) &
            +f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m9
    m(9) = 2.0d0*(f(1,i,j,k)+f(2,i,j,k))-f(3,i,j,k)-f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m10
     m(10) = -4.0d0*(f(1,i,j,k)+f(2,i,j,k))+2.0d0*(f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k)) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m11
    m(11) = f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m12
     m(12) = -2.0d0*(f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k))+( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k) )-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m13
    m(13) = f(7,i,j,k)-f(8,i,j,k)-f(9,i,j,k)+f(10,i,j,k)
    !---m14
    m(14) = f(15,i,j,k)-f(16,i,j,k)-f(17,i,j,k)+f(18,i,j,k)
    !---m15
    m(15) = f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m16
    m(16) = f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)-f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m17
    m(17) = -f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m18
    m(18) = f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)-f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !--------------------------------------------------------------------------------------------------------------------

    meq(0) = rho(i,j,k)
    meq(1) = -11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(2) = 3.0d0*rho(i,j,k)-11.0d0/2.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(3) = rho(i,j,k)*u(i,j,k)
    meq(4) = -2.0d0/3.0d0*meq(3)
    meq(5) = rho(i,j,k)*v(i,j,k)
    meq(6) = -2.0d0/3.0d0*meq(5)
    meq(7) = rho(i,j,k)*w(i,j,k)
    meq(8) = -2.0d0/3.0d0*meq(7)
    meq(9) = rho(i,j,k)*(2.0d0*u(i,j,k)*u(i,j,k)-v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(10) = -0.5d0*meq(9)
    meq(11) = rho(i,j,k)*(v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(12) = -0.5d0*meq(11)
    meq(13) = rho(i,j,k)*(u(i,j,k)*v(i,j,k))
    meq(14) = rho(i,j,k)*(v(i,j,k)*w(i,j,k))
    meq(15) = rho(i,j,k)*(w(i,j,k)*u(i,j,k))
    meq(16) = 0.0d0
    meq(17) = 0.0d0
    meq(18) = 0.0d0

            s(0) = 0.0d0 
            s(1) = Snu  !!!s_{e}
            s(2) = Snu   !!! s_{epsilon}
            s(3) = 0.0d0 
            s(4) = Sq   !!! s_{q}
            s(5) = 0.0d0 
            s(6) = Sq   !!! s_{q}
            s(7) = 0.0d0 
            s(8) = Sq   !!! s_{q}
            s(9) = Snu !!! s_{nu}
            s(10) = Snu   !!! s_{pi}
            s(11) = Snu   !!! s_{nu}
            s(12) = Snu !!! s_{pi}
            s(13) = Snu !!! s_{nu}
            s(14) = Snu   !!! s_{nu}
            s(15) = Snu   !!! s_{nu}
            s(16) = Sq   !!! s_{m}
            s(17) = Sq   !!! s_{m}
            s(18) = Sq   !!! s_{m}

            Fx(i,j,k) = -2.0d0*rho(i,j,k)*v(i,j,k)*omegaRatating
            Fy(i,j,k) = 2.0d0*rho(i,j,k)*u(i,j,k)*omegaRatating
            Fz(i,j,k) = rho(i,j,k)*gBeta*(T(i,j,k)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = 38.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(2) = -11.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(3) = Fx(i,j,k)
            fSource(4) = -2.0d0/3.0d0*Fx(i,j,k)
            fSource(5) = Fy(i,j,k)
            fSource(6) = -2.0d0/3.0d0*Fy(i,j,k)
            fSource(7) = Fz(i,j,k)
            fSource(8) = -2.0d0/3.0d0*Fz(i,j,k)
            fSource(9) = 4.0d0*u(i,j,k)*Fx(i,j,k)-2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(10) = -2.0d0*u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k) 
            fSource(11) = 2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(12) = -v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k)
            fSource(13) = u(i,j,k)*Fy(i,j,k)+v(i,j,k)*Fx(i,j,k)
            fSource(14) = v(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fy(i,j,k)
            fSource(15) = u(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fx(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+(1.0d0-0.5d0*s(alpha))*fSource(alpha)
            enddo

    f_post(0,i,j,k) = m_post(0)/19.0d0-5.0d0/399.0d0*m_post(1)+m_post(2)/21.0d0

    f_post(1,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(3,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(4,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(5,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

!---
    f_post(7,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0+( m_post(16)-m_post(17) )*0.125d0

    f_post(8,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0-( m_post(16)+m_post(17) )*0.125d0

    f_post(9,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0+( m_post(16)+m_post(17) )*0.125d0

    f_post(10,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0-( m_post(16)-m_post(17) )*0.125d0

!---
    f_post(11,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)-0.1250d0*( m_post(16)-m_post(18) )

    f_post(12,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)+0.125d0*( m_post(16)+m_post(18) )

    f_post(13,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)-0.125d0*( m_post(16)+m_post(18) )

    f_post(14,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)+0.125d0*( m_post(16)-m_post(18) )

!---
    f_post(15,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)+0.125d0*( m_post(17)-m_post(18) )

    f_post(16,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)-0.125d0*( m_post(17)+m_post(18) )

    f_post(17,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)+0.125d0*( m_post(17)+m_post(18) )

    f_post(18,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)-0.125d0*( m_post(17)-m_post(18) )

            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine collision


    subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,18
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine streaming


    subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef noslipWalls
    !$omp parallel do default(none) shared(f,f_post) private(j,k)
    do k=1,nz
        do j=1,ny
            !Back plane (i = 1)
            f(1,1,j,k) = f_post(2,1,j,k)
            f(7,1,j,k) = f_post(10,1,j,k)
            f(9,1,j,k) = f_post(8,1,j,k)
            f(11,1,j,k) = f_post(14,1,j,k)
            f(13,1,j,k) = f_post(12,1,j,k)

            !Front plane (i=nx)
            f(2,nx,j,k) = f_post(1,nx,j,k)
            f(8,nx,j,k) = f_post(9,nx,j,k)
            f(10,nx,j,k) = f_post(7,nx,j,k)
            f(12,nx,j,k) = f_post(13,nx,j,k)
            f(14,nx,j,k) = f_post(11,nx,j,k)
        enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do default(none) shared(f,f_post) private(i,k)
    do k=1,nz
        do i=1,nx
            !Left plane (j=1)
            f(3,i,1,k) = f_post(4,i,1,k)
            f(7,i,1,k) = f_post(10,i,1,k)
            f(8,i,1,k) = f_post(9,i,1,k)
            f(15,i,1,k) = f_post(18,i,1,k)
            f(17,i,1,k) = f_post(16,i,1,k)

            !Right plane (j=ny)
            f(4,i,ny,k) = f_post(3,i,ny,k)
            f(9,i,ny,k) = f_post(8,i,ny,k)
            f(10,i,ny,k) = f_post(7,i,ny,k)
            f(16,i,ny,k) = f_post(17,i,ny,k)
            f(18,i,ny,k) = f_post(15,i,ny,k)
        enddo
    enddo
    !$omp end parallel do
#endif

    !$omp parallel do default(none) shared(f,f_post) private(i,j)
    do j=1,ny
        do i=1,nx
            !Bottom side (k=1)
            f(5,i,j,1) = f_post(6,i,j,1)
            f(11,i,j,1) = f_post(14,i,j,1)
            f(12,i,j,1) = f_post(13,i,j,1)
            f(15,i,j,1) = f_post(18,i,j,1)
            f(16,i,j,1) = f_post(17,i,j,1)

            !Top side (k=nz)
            f(6,i,j,nz) = f_post(5,i,j,nz)
            f(13,i,j,nz) = f_post(12,i,j,nz)
            f(14,i,j,nz) = f_post(11,i,j,nz)
            f(17,i,j,nz) = f_post(16,i,j,nz)
            f(18,i,j,nz) = f_post(15,i,j,nz)
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine bounceback
    

    subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(f,rho,u,v,w,Fx,Fy,Fz) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k) = f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
            
                u(i,j,k) = ( f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)+0.5d0*Fx(i,j,k) )/rho(i,j,k)
            
                v(i,j,k) = ( f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fy(i,j,k) )/rho(i,j,k)
            
                w(i,j,k) = ( f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fz(i,j,k) )/rho(i,j,k)
            
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine macro


    subroutine collisionT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !------------------------
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: q(0:6)

    !$omp parallel do default(none) shared(g,g_post,u,v,w,T) private(i,j,k,alpha,n,neq,q,n_post) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
            
    n(0) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(1) = g(1,i,j,k)-g(2,i,j,k)
    n(2) = g(3,i,j,k)-g(4,i,j,k)
    n(3) = g(5,i,j,k)-g(6,i,j,k)
    n(4) = -6.0d0*g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(5) = 2.0d0*g(1,i,j,k)+2.0d0*g(2,i,j,k)-g(3,i,j,k)-g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
    n(6) = g(3,i,j,k)+g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
        
            neq(0) = T(i,j,k)
            neq(1) = T(i,j,k)*u(i,j,k)
            neq(2) = T(i,j,k)*v(i,j,k)
            neq(3) = T(i,j,k)*w(i,j,k)
            neq(4) = T(i,j,k)*paraA
            neq(5) = 0.0d0
            neq(6) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qd
            q(4) = Qnu
            q(5) = Qnu
            q(6) = Qnu
        
            do alpha=0,6
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j,k) = n_post(0)/7.0d0-n_post(4)/7.0d0
    g_post(1,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(2,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(3,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(4,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(5,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
    g_post(6,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
        
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine collisionT


    subroutine streamingT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(g,g_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,6
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    g(alpha,i,j,k) = g_post(alpha,ip,jp,kp)
                    
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine streamingT


    subroutine bouncebackT()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef TopBottomPlatesAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(i,j)
    do j=1,ny
        do i=1,nx
            !Bottom side (k=1)
            g(5,i,j,1) = g_post(6,i,j,1)

            !Top side (k=nz)
            g(6,i,j,nz) = g_post(5,i,j,nz)
        enddo
    enddo
    !$omp end parallel do
#endif

#ifdef TopBottomPlatesConstT
    !$omp parallel do default(none) shared(g,g_post) private(i,j)
    do j=1,ny
        do i=1,nx
            !Bottom side (k=1)
            g(5,i,j,1) = -g_post(6,i,j,1)+(6.0d0+paraA)/21.0d0*Thot

            !Top side (k=nz)
            g(6,i,j,nz) = -g_post(5,i,j,nz)+(6.0d0+paraA)/21.0d0*Tcold
        enddo
    enddo
    !$omp end parallel do
#endif


#ifdef LeftRightWallsConstT
    !$omp parallel do default(none) shared(g,g_post) private(i,k)
    do k=1,nz
        do i=1,nx
            !Left side (j=1)
            g(3,i,1,k) = -g_post(4,i,1,k)+(6.0d0+paraA)/21.0d0*Thot

            !Right side (j=ny)
            g(4,i,ny,k) = -g_post(3,i,ny,k)+(6.0d0+paraA)/21.0d0*Tcold
        enddo
    enddo
    !$omp end parallel do
#endif

#ifdef LeftRightWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(i,k)
    do k=1,nz
        do i=1,nx
            !Left side (j=1)
            g(3,i,1,k) = g_post(4,i,1,k)

            !Right side (j=ny)
            g(4,i,ny,k) = g_post(3,i,ny,k)
        enddo
    enddo
    !$omp end parallel do
#endif
    
    
#ifdef BackFrontWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(j,k)
    do k=1,nz
        do j=1,ny
            !Front side (i=nx)
            g(2,nx,j,k) = g_post(1,nx,j,k)

            !Back side (i=1)
            g(1,1,j,k) = g_post(2,1,j,k)
        enddo
    enddo
    !$omp end parallel do
#endif

    return
    end subroutine bouncebackT


    subroutine macroT()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(g,T) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                T(i,j,k) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine macroT


    subroutine check()
    use commondata
    implicit none
    integer :: i, j, k
    real(kind=8) :: error1, error2, error5, error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,w,wp,T,Tp) private(i,j,k) reduction(+:error1,error2,error5,error6)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                error1 = error1+(u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k))+(v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k))+(w(i,j,k)-wp(i,j,k))*(w(i,j,k)-wp(i,j,k))
                error2 = error2+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                
                error5 = error5+dABS( T(i,j,k)-Tp(i,j,k) )
                error6 = error6+dABS( T(i,j,k) )
                
                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
                wp(i,j,k) = w(i,j,k)
                Tp(i,j,k) = T(i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do 
    
    errorU = dsqrt(error1)/dsqrt(error2)
    errorT = error5/error6

    write(*,*) itc,' ',errorU,' ',errorT

    return
    end subroutine check
    
    
    subroutine output_binary()
    use commondata
    implicit none
    integer :: i, j, k
    character(len=100) :: filename

#ifdef RBconvection
    fileNum = fileNum+1
    write(filename,*) fileNum
    filename = adjustl(filename)
#endif
#ifdef benchmarkCavity
    write(filename,*) itc
    filename = adjustl(filename)
#endif

    open(unit=01,file='buoyancyCavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(01)

    return
    end subroutine output_binary
    
    
    subroutine backupData()
    use commondata
    implicit none
    integer :: i, j, k, alpha
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=01,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) ((((f(alpha,i,j,k),alpha=0,18),i=1,nx),j=1,ny),k=1,nz)
    write(01) ((((g(alpha,i,j,k),alpha=0,6),i=1,nx),j=1,ny),k=1,nz)
    close(01)

    return
    end subroutine backupData


#ifdef outputTecplot
    subroutine output_Tecplot()
    use commondata
    implicit none
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    character(len=100) :: filename
 
#ifdef RBconvection
    fileNum = fileNum+1
    write(filename,*) fileNum
    filename = adjustl(filename)
#endif
#ifdef benchmarkCavity
    write(filename,*) itc
    filename = adjustl(filename)
#endif
    open(41,file='buoyancyCavity-'//trim(filename)//'.plt', access='stream', form='unformatted')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--Integer value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file (here 5 variables)
    write(41) 7

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='Z'
    call dumpstring(V3)
    
    V4='U'
    call dumpstring(V4)
    V5='V'
    call dumpstring(V5)
    V6='W'
    call dumpstring(V6)
    
    V7='T'
    call dumpstring(V7)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName='ZONE 001'
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
                write(41) real(T(i,j,k))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
    end subroutine output_Tecplot


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
#endif

    subroutine calNuRe()
    use commondata
    implicit none
    integer :: i, j, k 
    integer :: halfTime
    real(kind=8) :: NuVolAvg_firstHalf, NuVolAvg_lastHalf
    real(kind=8) :: ReVolAvg_firstHalf, ReVolAvg_lastHalf
    real(kind=8) :: NuVolAvg_temp
    real(kind=8) :: ReVolAvg_temp
    
    dimensionlessTime = dimensionlessTime+1
    !-------------------------------------------------------------------------------------------------------
    NuVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(w,T) private(i,j,k) reduction(+:NuVolAvg_temp)
    do k=1,nz
        do j=1,ny
            do i=1,nx       
                NuVolAvg_temp = NuVolAvg_temp+w(i,j,k)*T(i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(nx*ny*nz)*dble(nz)/diffusivity+1.0d0
    
    open(unit=01,file="Nu_Bulk.dat",status='unknown',position='append')
    write(01,*) itc/timeUnit, NuVolAvg(dimensionlessTime)
    close(01)
    
    NuVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)+NuVolAvg(i) 
    enddo
    NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------
    ReVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(u,v,w) private(i,j,k) reduction(+:ReVolAvg_temp)
    do k=1,nz
        do j=1,ny
            do i=1,nx       
                ReVolAvg_temp = ReVolAvg_temp+(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))
            enddo
        enddo
    enddo
    !$omp end parallel do
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(nx*ny*nz))*dble(nz)/viscosity
    
    open(unit=02,file="./Re_Bulk.dat",status='unknown',position='append')
    write(02,*) itc/timeUnit, ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)

    ReVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)+ReVolAvg(i) 
    enddo
    ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------    

    if(statisticallyStationaryState.EQ.1) then
        if(MOD(dimensionlessTime-statisticallyStationaryTime,200).EQ.0) then
            
            halfTime = (statisticallyStationaryTime+dimensionlessTime)/2
            
            NuVolAvg_firstHalf = 0.0d0
            ReVolAvg_firstHalf = 0.0d0
            do i=statisticallyStationaryTime,halfTime
                NuVolAvg_firstHalf = NuVolAvg_firstHalf+NuVolAvg(i)
                ReVolAvg_firstHalf = ReVolAvg_firstHalf+ReVolAvg(i)
            enddo
            NuVolAvg_firstHalf = NuVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
            ReVolAvg_firstHalf = ReVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
            
            NuVolAvg_lastHalf = 0.0d0
            ReVolAvg_lastHalf = 0.0d0
            do i=halfTime,dimensionlessTime
                NuVolAvg_lastHalf = NuVolAvg_lastHalf+NuVolAvg(i)
                ReVolAvg_lastHalf = ReVolAvg_lastHalf+ReVolAvg(i)
            enddo
            NuVolAvg_lastHalf = NuVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
            ReVolAvg_lastHalf = ReVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
            
            if( (dabs((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf).LE.0.01d0).AND. &
                 (dabs((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf).LE.0.01d0) ) then

                if(flowReversal.EQ.0) then
                    statisticallyConvergeState = 1
                elseif(flowReversal.EQ.1) then
                    if(dimensionlessTime.GE.(statisticallyStationaryTime+flowReversalTime)) statisticallyConvergeState = 1
                endif
                
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|---------Statistically Converge Reached at Time =",dimensionlessTime,"----------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif
            
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "statisticallyStationaryTime=",statisticallyStationaryTime
            write(01,*) "halfTime =",halfTime
            write(01,*) "dimensionlessTime =",dimensionlessTime
            write(01,*) "NuVolAvg diff. (first and last half):",dabs(((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (first and last half):",dabs(((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "               "
            close(01)        
        endif
    endif

    if(statisticallyStationaryState.EQ.0) then
        if(MOD(dimensionlessTime,100).EQ.0) then
            if( (dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-100))/NuVolAvg_mean(dimensionlessTime)).LE.0.01d0).AND. &
                (dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-100))/ReVolAvg_mean(dimensionlessTime)).LE.0.01d0) ) then
                statisticallyStationaryState = 1
                statisticallyStationaryTime = dimensionlessTime
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|-----------Statistically Stationary Reached at Time =",statisticallyStationaryTime,"------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                write(01,*) "   "
                close(01)
            endif 
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =",dimensionlessTime
            write(01,*) "NuVolAvg diff. (100 timeunit interval):",dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-100)) &
                                                                            /NuVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (100 timeunit interval):",dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-100)) &
                                                                            /ReVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "               "
            close(01)       
        endif
    endif
    
    return
    end subroutine calNuRe