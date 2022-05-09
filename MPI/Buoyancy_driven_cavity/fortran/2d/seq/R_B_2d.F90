!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


!~!~For steadyFlow, we define array up, vp, Tp to check  convergence;
!~!~for unsteadyFlow, we did not define array up, vp, Tp to save memory.
! #define steadyFlow    
#define unsteadyFlow

!~!~Uncomment below to simulate mass particles
!~#define pointParticle

!!!~~velocity B.C.~~
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!!!~#define VerticalWallsPeriodicalU
!!!~#define HorizontalWallsFreeslip
!!!~#define VerticalWallsFreeslip
!!!~~velocity B.C.~~

!!!!~~temperature B.C. (for Rayleigh Benard Cell)~~
#define RayleighBenardCell
#define HorizontalWallsConstT
#define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT
!!!!~~temperature B.C.~~

!!!~~temperature B.C. (for Side Heated Cell)~~
! #define SideHeatedCell
! #define HorizontalWallsAdiabatic
! #define VerticalWallsConstT
!!!~~temperature B.C.~~


    module ioFolder
        character(len=100) :: binFolderPrefix="../binFile/buoyancyCavity"
        character(len=100) :: pltFolderPrefix="../pltFile/buoyancyCavity"
        character(len=100) :: porousGeoFile="../case1/obstSquareCylinder-126.bin"
        character(len=100) :: particlePositionFolderPrefix="../particlePositionFile/particlePosition"
        character(len=100) :: particleReynoldsFolderPrefix="../particleReynoldsFile/particleReynolds"
    end module ioFolder
    
    module commondata
        implicit none
        !-----------------------------------------------------------------------------------------------
        integer(kind=4), parameter :: loadInitField=0 !! value 0: do NOT load previous data file; value 1: load data from previous one.
        integer(kind=4), parameter :: reloadDimensionlessTime=0 ! total free-fall time in previous simulation, check Re_Vol.dat // Nu_Bulk.dat
        integer(kind=4), parameter :: reloadStatisticallyConvergeTime=0 !free-fall time after reaching StatisticallyConvergeState in previous simulation.
        integer(kind=4), parameter :: reloadbinFileNum=0 ! check backup-*.bin file name 
        
        !! value 0: do NOT simulate reversal; value 1: do simulate reversal
        !! for reversal, the simulation stops at "t = flowReversalTimeMax" (after stationary)
        !! for non-reversal, the simulation stops when statistics converge.
        integer(kind=4), parameter :: flowReversal=0  
        
        integer(kind=4), parameter :: nx=201, ny=201     !----Section 1----!
        integer(kind=4), parameter :: flowGeometry=1    !----Section 1----!
        real(kind=8), parameter :: lengthUnit=dble(ny)     !----Section 1----!
        
        real(kind=8), parameter :: Rayleigh=1e7        !----Section 2----!
        real(kind=8), parameter :: Prandtl=5.3d0       !----Section 2----!
        real(kind=8), parameter :: Mach=0.1d0           !----Section 2----!
        
        real(kind=8), parameter :: outputFrequency=1.0d0 !~unit free fall time                            !----Section 3----!
        
        integer(kind=4), parameter :: dimensionlessTimeMax=int(3000/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: flowReversalTimeMax=int(10000/outputFrequency)    !----Section 3----!
                                                        !~!~if flowReversal=1,dimensionlessTimeMax shoud > flowReversalTimeMax
        integer(kind=4), parameter :: backupInterval=2000        !~!~unit: free-fall time            !----Section 3----!
        integer(kind=4), parameter :: minAvgPeriod=int(1000/outputFrequency)               !----Section 3----!
        integer(kind=4), parameter :: stationaryCheckInterval=int(100/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: convergeCheckInterval=int(100/outputFrequency)  !----Section 3----!
        
        real(kind=8), parameter :: statStationError=0.02d0  !----Section 3----!
        real(kind=8), parameter :: statConvError=0.02d0     !----Section 3----!
        real(kind=8), parameter :: epsU=1e-6                     !----Section 3----!
        real(kind=8), parameter :: epsT=1e-6                     !----Section 3----!

        real(kind=8), parameter :: shearReynolds=100.0d0                                                          !----Section 4----!
        
        integer(kind=4), parameter :: outputBinFile=0, outputPltFile=1                                              !----Section *----!
                
        real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
        real(kind=8), parameter :: tiltedAngle=dble(0.0d0/180.0d0*Pi)
        !-----------------------------------------------------------------------------------------------        
        !----Section 1----!
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        integer(kind=4), parameter :: nxFourth=(nx-1)/4+1, nyFourth=(ny-1)/4+1
        integer(kind=4), parameter :: nxThreeFourths=3*(nx-1)/4+1, nyThreeFourths=3*(ny-1)/4+1
        real(kind=8), parameter :: xCenter=dble(nxHalf), yCenter=dble(nyHalf)
        
        integer(kind=4), parameter :: squareGeo=1, cornerLessGeo=2, circularGeo=3, porousIsoGeo=4
       
        !-----------------------------------------------------------------------------------------------
        !----Section 2----!
        real(kind=8), parameter :: rho0=1.0d0 !~!~m.u.
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8), parameter :: tauf=0.5d0+Mach*lengthUnit*DSQRT(3.0d0*Prandtl/Rayleigh)
        real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
        real(kind=8), parameter :: diffusivity=viscosity/Prandtl
        
        real(kind=8), parameter :: paraA=20.0d0*dsqrt(3.0d0)*diffusivity-4.0d0
        real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/lengthUnit
        real(kind=8), parameter :: gBeta=gBeta1/lengthUnit/lengthUnit
        
        real(kind=8), parameter :: timeUnit=dsqrt(lengthUnit/gBeta)  !!dble(ny*ny)/diffusivity
        real(kind=8), parameter :: velocityUnit=dsqrt(gBeta*lengthUnit)
    
        real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
        !-----------------------------------------------------------------------------------------------
        !----Section 3----!
        integer(kind=4) :: statisticallyStationaryState, statisticallyStationaryTime
        integer(kind=4) :: statisticallyConvergeState
        real(kind=8) :: errorU, errorT
        !-----------------------------------------------------------------------------------------------
        !----Section 4----!
        real(kind=8), parameter :: U0=shearReynolds*viscosity/dble(ny)
        real(kind=8), parameter :: UwallTopLeft=U0, UwallTopRight=-U0, UwallBottomLeft=-U0, UwallBottomRight=U0
        real(kind=8), parameter :: UwallLeftTop=U0, UwallLeftBottom=U0, UwallRightTop=U0, UwallRightBottom=U0  !----Section 4----!!----Section 4----!
        !-----------------------------------------------------------------------------------------------
        !----Section 5----!
        real(kind=8) :: xp(0:nx+1), yp(0:ny+1)
        real(kind=8), allocatable :: u(:,:), v(:,:), T(:,:), rho(:,:)
#ifdef steadyFlow
        real(kind=8), allocatable :: up(:,:), vp(:,:), Tp(:,:)
#endif
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)
        real(kind=8), allocatable :: g(:,:,:), g_post(:,:,:)
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)
        
        integer(kind=4) :: ex(0:8), ey(0:8)
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
        integer(kind=4) :: r(1:8)
        data r/3, 4, 1, 2, 7, 8, 5, 6/
        real(kind=8) :: omega(0:8), omegaT(0:4)
        !-----------------------------------------------------------------------------------------------
        !----Section 6----!
        integer(kind=4) :: itc
        integer(kind=4), parameter :: itc_max=dimensionlessTimeMax*int(outputFrequency*timeUnit)
        
        integer(kind=4) :: binFileNum, pltFileNum
        integer(kind=4) :: particleFileNum
        integer(kind=4) :: dimensionlessTime
        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)  !*here, the index should start from zero
        real(kind=8) :: NuVolAvg_mean(0:dimensionlessTimeMax), ReVolAvg_mean(0:dimensionlessTimeMax)  !*because we will calculate Nu(t)-Nu(t-100) at time t=100.
        !-----------------------------------------------------------------------------------------------
        !----Section 7----!
        integer(kind=4) :: fluidNumMax
        integer(kind=4), allocatable :: obst(:,:)
        !-----------------------------------------------------------------------------------------------        
    end module commondata
        
    
!*!*Begin program main
    program main
    use omp_lib    
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    real(kind=8) :: timeStart2, timeEnd2
    integer(kind=4) :: myMaxThreads
    INTEGER(kind=4) :: time
    character(len=24) :: ctime, string
        
    open(unit=00,file="SimulationSettings.txt",status='unknown',position="append")
    string = ctime( time() )
    write(00,*) 'Start: ', string
#ifdef _OPENMP
    write(00,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(24)
    myMaxThreads = OMP_get_max_threads()
    write(00,*) "Max Running threads:",myMaxThreads
#endif
    close(00)
    
    call initial()

    if(outputPltFile.EQ.1) then
        call output_Tecplot()
    endif
                
    call CPU_TIME(timeStart)
#ifdef _OPENMP
    timeStart2 = OMP_get_wtime()
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

#ifdef steadyFlow
        if(MOD(itc,2000).EQ.0) call check()
#endif

#ifdef SideHeatedCell
        if( (MOD(itc,backupInterval*int(timeUnit)).EQ.0).AND.(outputBinFile.EQ.1) ) call backupData()
#endif


#ifdef RayleighBenardCell
        if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then
            
#ifdef steadyFlow
            if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()
#endif
            
            if(statisticallyStationaryState.EQ.1) then
                if(outputBinFile.EQ.1) then
                    call output_binary()
                    if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()
                endif   
                if(outputPltFile.EQ.1) then
                    call output_Tecplot()
                endif
            endif
            
            call calNuRe()
            
        endif
#endif

    enddo
    
    if(outputBinFile.EQ.1) call backupData()
    
    if(outputPltFile.EQ.1) call output_Tecplot()
        
    call CPU_TIME(timeEnd)
#ifdef _OPENMP
    timeEnd2 = OMP_get_wtime()
#endif
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
#ifdef _OPENMP
    write(00,*) "Time (OMP) = ", real(timeEnd2-timeStart2), "s"
#endif
    
    write(00,*) "Dellocate Array......"
    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(T)
#ifdef steadyFlow
    deallocate(up)
    deallocate(vp)
    deallocate(Tp)
#endif
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(obst)
    write(00,*) "    "

    write(00,*) "Successfully: DNS completed!"
    
    string = ctime( time() )
    write(00,*) 'End:   ', string
    close(00)
    
    stop
    end program main
!*!*end program main


!*!*Begin subroutine initial    
    subroutine initial()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2
    character(len=100) :: reloadFileName
    real(kind=8) :: dev
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    fluidNumMax = nx*ny
    
open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')

    if(outputBinFile.EQ.1) then
        open(unit=01,file=trim(binFolderPrefix)//"-"//"readme",status="unknown")
        write(01,*) "binFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", binFolderPrefix
    endif 
    if(outputPltFile.EQ.1) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    
    if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then
        write(00,*) '----------------------------------'
        write(00,*) 'Error: paraA=', paraA
        write(00,*) '... Current Mach number is :', real(Mach)
        write(00,*) '... Please try to reduce Mach number!'
        write(00,*) '----------------------------------'
        stop
    endif

    if(flowReversal.EQ.1) then
        if(dimensionlessTimeMax.LE.flowReversalTimeMax) then
            write(00,*) "Error: in case of flor reversal,"
            write(00,*) "... please check dimensionlessTimeMax!"
            stop
        endif
    endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',nx,ny
    write(00,*) 'Rayleigh=',real(Rayleigh), '; Prandtl =',real(Prandtl), '; Mach =',real(Mach)
    write(00,*) 'shearReynolds=', real(shearReynolds)
    if(shearReynolds.GT.0.0d0) then
        write(00,*) "Richardson=", real(Rayleigh/(shearReynolds**2.0d0*Prandtl))
        write(00,*) "...UwallTopLeft=", real(UwallTopLeft), "; UwallTopRight=", real(UwallTopRight)
        write(00,*) "...UwallBottomLeft=", real(UwallBottomLeft), "; UwallBottomRight=", real(UwallBottomRight)
        write(00,*) "...UwallLeftTop=", real(UwallLeftTop), "; UwallLeftBottom=", real(UwallLeftBottom)
        write(00,*) "...UwallRightTop=", real(UwallRightTop), "; UwallRightBottom=", real(UwallRightBottom)
    endif
    write(00,*) "Length unit: L0 =", real(lengthUnit)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf)
    write(00,*) "viscosity =",real(viscosity), "; diffusivity =",real(diffusivity)
    write(00,*) "outputFrequency", real(outputFrequency), "tf"
    write(00,*) "......................  or ",  int(outputFrequency*timeUnit), "in itc units"
    write(00,*) "dimensionlessTimeMax:"
    write(00,*) "...statistically stationary PLUS convergence = ", real(dimensionlessTimeMax*outputFrequency)
    if(flowReversal.EQ.1) then
        write(00,*) "flowReversalTimeMax:"
        write(00,*) "...statistically stationary PLUS convergence = ", real(flowReversalTimeMax*outputFrequency)
    endif
    write(00,*) "minAvgPeriod:"
    write(00,*) "...to obtain statistically convergence =", real(minAvgPeriod*outputFrequency)
    write(00,*) "backupInterval =", backupInterval, " free-fall time units"
    write(00,*) ".................... or ", int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit), "itc units"
    if(loadInitField.EQ.1) then
        write(00,*) "reloadDimensionlessTime=", reloadDimensionlessTime
        write(00,*) "reloadStatisticallyConvergeTime=", reloadStatisticallyConvergeTime
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU),"; epsT =", real(epsT)
    write(00,*) "    "
    

    write(00,*) "I am regular geometry (Square Cavity)"
    if(flowGeometry.NE.squareGeo) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "... flowGeometry should be", squareGeo, "; while its actual value is ", flowGeometry
        stop
    endif
    if(dabs(lengthUnit-dble(ny)).GT.1e-6) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "... lengthUnit should be equal to ny"
        stop
    endif

#ifdef RayleighBenardCell
    write(00,*) "I am Rayleigh Benard Cell"
#endif
#ifdef SideHeatedCell
    write(00,*) "I am Side Heated Cell"
#endif

#ifdef steadyFlow
        write(00,*) "I am steadyFlow"
#endif
#ifdef unsteadyFlow
        write(00,*) "I am unsteadyFlow"
#endif
    
    if(flowReversal.EQ.1) then
        write(00,*) "I am also flow reversal"
    elseif(flowReversal.EQ.0) then
        write(00,*) "I am NOT flow reversal"
    else
        write(00,*) "Error: Please check flowReversal setting!"
        stop
    endif

    write(00,*) "I am level cell"


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
    allocate (T(nx,ny))
    allocate (rho(nx,ny))
    
#ifdef steadyFlow
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (Tp(nx,ny))
#endif
    
    allocate (f(0:8,nx,ny))
    allocate (f_post(0:8,0:nx+1,0:ny+1))
    allocate (g(0:4,nx,ny))
    allocate (g_post(0:4,0:nx+1,0:ny+1))
    
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))
    
    allocate (obst(0:nx+1,0:ny+1))

    obst = 0
    rho = rho0
    
    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/5.0d0
    do alpha=1,4
        omegaT(alpha) = (paraA+4.0d0)/20.0d0
    enddo
    

    if(loadInitField.EQ.0) then 
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        
        do i=1,nxHalf
            u(i,ny) = UwallTopLeft
            u(i,1)  = UwallBottomLeft
        enddo
        do i=nxHalf+1, nx
            u(i,ny) = UwallTopRight
            u(i,1) = UwallBottomRight
        enddo
        do j=2,nyHalf
            v(1,j) = UwallLeftBottom
            v(nx,j) = UwallRightBottom
        enddo
        do j=nyHalf+1, ny-1
            v(1,j) = UwallLeftTop
            v(nx,j) = UwallRightTop
        enddo
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0) then
            write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
            stop
        endif
        
#ifdef VerticalWallsNoslip
        write(00,*) "Velocity B.C. for vertical walls are: ===No-slip wall==="
#endif
#ifdef VerticalWallsFreeslip
        write(00,*) "Velocity B.C. for vertical walls are: ===Free-slip wall==="
#endif
#ifdef VerticalWallsPeriodicalU
        write(00,*) "Velocity B.C. for vertical walls are: ===Periodical==="
#endif
#ifdef HorizontalWallsNoslip
        write(00,*) "Velocity B.C. for horizontal walls are: ===No-slip wall==="
#endif
#ifdef HorizontalWallsFreeslip
        write(00,*) "Velocity B.C. for horizontal walls are: ===Free-slip wall==="
#endif



#ifdef VerticalWallsConstT
    do j=1,ny
        !~T(1,j) = Thot
        !~T(nx,j) = Tcold
        do i=1,nx
            T(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
        enddo
    enddo
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif
#ifdef HorizontalWallsConstT
    do i=1,nx
        !~T(i,1) = Thot
        !~T(i,ny) = Tcold
        do j=1,ny
            T(i,j) = dble(j-1)/dble(ny-1)*(Tcold-Thot)+Thot
        enddo
    enddo
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
#endif

#endif
    
#ifdef VerticalWallsPeriodicalT
    write(00,*) "Temperature B.C. for vertical walls are:===Periodical wall==="
#endif

#ifdef VerticalWallsAdiabatic
    write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif
        
        f = 0.0d0
        g = 0.0d0
        !$omp parallel do default(none) shared(f,g,u,v,T,ex,ey,omega,omegaT,rho) private(i,j,alpha,us2,un)
        do j=1,ny
            do i=1,nx
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha=0,8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo
                do alpha=0,4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    g(alpha,i,j) = T(i,j)*omegaT(alpha)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))
                enddo
            enddo
        enddo
        !$omp end parallel do

    elseif(loadInitField.EQ.1) then
        if(reloadDimensionlessTime.EQ.0) then
            write(00,*) "Error: since loadInitField.EQ.1, reloadDimensionlessTime should not be 0"
            stop
        endif
        write(00,*) "Load initial field from previous simulation >>>"
        write(reloadFileName, *) reloadbinFileNum
        reloadFileName = adjustl(reloadFileName)
        open(unit=01,file="../reloadFile/backupFile-"//trim(reloadFileName)//".bin",form="unformatted",access="sequential",status='old')
        read(01) (((f(alpha,i,j), alpha=0,8), i=1,nx), j=1,ny)
        read(01) (((g(alpha,i,j), alpha=0,4), i=1,nx), j=1,ny)
        read(01) ((u(i,j),i=1,nx), j=1,ny)
        read(01) ((v(i,j),i=1,nx), j=1,ny)
        reaD(01) ((T(i,j),i=1,nx), j=1,ny)
        close(01)
        
        dev = 0.0d0
        do j=1,ny
            do i=1,nx
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                dev = dev+g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            enddo
        enddo
        write(00,*) "RELOAD: Deviation in temperature: ", real(dev)
        if(dev.GT.1.0d0) then
            write(00,*) "Error: too large Deviation when reload data!"
            stop
        endif
        write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"
    else
        write(00,*) "Error: initial field is not properly set"
        stop
    endif
    
    write(00,*)"-------------------------------------------------------------------------------"
close(00)

#ifdef steadyFlow
    up = 0.0d0
    vp = 0.0d0
    Tp = 0.0d0
#endif
        
    f_post = 0.0d0
    g_post = 0.0d0
    
    binFileNum = 0
    pltFileNum = 0
    particleFileNum = 0
    dimensionlessTime = 0

    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 1
    statisticallyStationaryTime = 0
#ifdef unsteadyFlow
    if(loadInitField.EQ.1) statisticallyStationaryState = 1
#endif
    statisticallyConvergeState = 0
    
    return
    end subroutine initial
!*!*end subroutine initial 
    
    
!*!*Begin subroutine collision
    subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,Fx,Fy,T,obst,omega) private(i,j,alpha,s,m,m_post,meq,fSource) 
    do j = 1, ny
        do i = 1, nx


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

            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))

            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0


        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine collision
!*!*End subroutine collision


!*!*Begin subroutine streaming
    subroutine streaming()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey,obst) private(i,j,ip,jp,alpha)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=0,8
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                
                    f(alpha,i,j) = f_post(alpha,ip,jp)
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine streaming
!*!*End subroutine streaming



!*!*Begin subroutine bounceback
    subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0, y0, q
    

#ifdef VerticalWallsPeriodicalU
    !$omp parallel do default(none) shared(f, f_post) private(j)
    do j=1,ny
        !Left side (i=1)
        f(1,1,j) = f_post(1,nx,j)
        f(5,1,j) = f_post(5,nx,j)
        f(8,1,j) = f_post(8,nx,j)
        
        !Right side (i=nx)
        f(3,nx,j) = f_post(3,1,j)
        f(6,nx,j) = f_post(6,1,j)
        f(7,nx,j) = f_post(7,1,j)
    enddo
    !$omp end parallel do
    i = 1
    j = 1
    f(2, i, j) = f_post(2, nx, j)
    f(6, i, j) = f_post(6, nx, j)
    
    i = nx
    j = 1
    f(2, i, j) = f_post(2, 1, j)
    f(5, i, j) = f_post(5, 1, j)
    
    i = 1
    j = ny
    f(4, i, j) = f_post(4, nx, j)
    f(7, i, j) = f_post(7, nx, j)
    
    i = nx
    j = ny
    f(4, i, j) = f_post(4, 1, j)
    f(8, i, j) = f_post(8, 1, j)
#endif

#ifdef VerticalWallsNoslip
    !$omp parallel do default(none) shared(f, f_post, rho) private(j)
    do j = 2, nyHalf
        !Left side (i=1)
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)-rho(1,j)*(-UwallLeftBottom)/6.0d0
        f(8,1,j) = f_post(6,1,j)-rho(1,j)*UwallLeftBottom/6.0d0
        
        !Right side (i=nx)
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)-rho(nx,j)*(-UwallRightBottom)/6.0d0
        f(7,nx,j) = f_post(5,nx,j)-rho(nx,j)*UwallRightBottom/6.0d0
    enddo
    !$omp end parallel do
    !$omp parallel do default(none) shared(f, f_post, rho) private(j)
    do j = nyHalf+1, ny-1
        !Left side (i=1)
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)-rho(1,j)*(-UwallLeftTop)/6.0d0
        f(8,1,j) = f_post(6,1,j)-rho(1,j)*UwallLeftTop/6.0d0
        
        !Right side (i=nx)
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)-rho(nx,j)*(-UwallRightTop)/6.0d0
        f(7,nx,j) = f_post(5,nx,j)-rho(nx,j)*UwallRightTop/6.0d0
    enddo
    !$omp end parallel do
#endif

!~ #ifdef VerticalWallsFreeslip
    !~!$omp parallel do default(none) shared(f, f_post) private(j)
    !~ do j=1,ny
        !~ !Left side (i=1)
        !~ f(1,1,j) = f(1,2,j)
        !~ f(5,1,j) = f(5,2,j)
        !~ f(8,1,j) = f(8,2,j)
        
        !~ !Right side (i=nx)
        !~ f(3,nx,j) = f(3,nx-1,j)
        !~ f(6,nx,j) = f(6,nx-1,j)
        !~ f(7,nx,j) = f(7,nx-1,j)
    !~ enddo
    !~!$omp end parallel do
!~ #endif

#ifdef HorizontalWallsNoslip
    !$omp parallel do default(none) shared(f, f_post, rho) private(i)
    do i = 2, nxHalf
        !Bottom side (j=1)
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)-rho(i,1)*(-UwallBottomLeft)/6.0d0
        f(6,i,1) = f_post(8,i,1)-rho(i,1)*UwallBottomLeft/6.0d0

        !Top side (j=ny)
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*UwallTopLeft/6.0d0
        f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-UwallTopLeft)/6.0d0
    enddo
    !$omp end parallel do
    !$omp parallel do default(none) shared(f, f_post, rho) private(i)
    do i = nxHalf+1, nx-1
        !Bottom side (j=1)
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)-rho(i,1)*(-UwallBottomRight)/6.0d0
        f(6,i,1) = f_post(8,i,1)-rho(i,1)*UwallBottomRight/6.0d0

        !Top side (j=ny)
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*UwallTopRight/6.0d0
        f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-UwallTopRight)/6.0d0
    enddo
    !$omp end parallel do
#endif

#ifndef VerticalWallsPeriodicalU
    !!~corners
    i = 1
    j = 1
    f(1,i,j) = f_post(3,i,j)
    f(2,i,j) = f_post(4,i,j)
    f(6,i,j) = f_post(8,i,j)-rho(i,j)*UwallBottomLeft/6.0d0
    f(8,i,j) = f_post(6,i,j)-rho(i,j)*UwallLeftBottom/6.0d0
    f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallBottomLeft-UwallLeftBottom)/6.0d0
    
    i = nx 
    j = 1
    f(2,i,j) = f_post(4,i,j)
    f(3,i,j) = f_post(1,i,j)
    f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallBottomRight)/6.0d0
    f(7,i,j) = f_post(5,i,j)-rho(i,j)*UwallRightBottom/6.0d0
    f(6,i,j) = f_post(8,i,j)-rho(i,j)*(UwallBottomRight-UwallRightBottom)/6.0d0

    i = 1
    j = ny 
    f(1,i,j) = f_post(3,i,j)
    f(4,i,j) = f_post(2,i,j)
    f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallLeftTop)/6.0d0
    f(7,i,j) = f_post(5,i,j)-rho(i,j)*UwallTopLeft/6.0d0
    f(8,i,j) = f_post(6,i,j)-rho(i,j)*(-UwallTopLeft+UwallLeftTop)/6.0d0
    
    i = nx 
    j = ny 
    f(3,i,j) = f_post(1,i,j)
    f(4,i,j) = f_post(2,i,j)
    f(6,i,j) = f_post(8,i,j)-rho(i,j)*(-UwallRightTop)/6.0d0
    f(8,i,j) = f_post(6,i,j)-rho(i,j)*(-UwallTopRight)/6.0d0
    f(7,i,j) = f_post(5,i,j)-rho(i,j)*(UwallTopRight+UwallRightTop)/6.0d0
    !!~corners
#endif

!~ #ifdef HorizontalWallsFreeslip
    !~!$omp parallel do default(none) shared(f) private(i)
    !~ do i=1,nx
        !~ !Bottom side (j=1)
        !~ f(2,i,1) = f(2,i,2)  
        !~ f(5,i,1) = f(5,i,2)  
        !~ f(6,i,1) = f(6,i,2)  

        !~ !Top side (j=ny)
        !~ f(4,i,ny) = f(4,i,ny-1)  
        !~ f(7,i,ny) = f(7,i,ny-1)  
        !~ f(8,i,ny) = f(8,i,ny-1)
    !~ enddo
    !~!$omp end parallel do
!~ #endif

    return
    end subroutine bounceback
!*!*End subroutine bounceback
    

!*!*Begin subroutine macro
    subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy,obst) private(i,j)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*Fx(i,j) )/rho(i,j)
                v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*Fy(i,j) )/rho(i,j)
            elseif(obst(i,j).EQ.1) then
                rho(i,j) = 1.0d0
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0
            endif
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine macro
!*!*End subroutine macro


!*!*Begin subroutine collisionT
    subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    !$omp parallel do default(none) shared(g,g_post,u,v,T,obst,omegaT) private(i,j,alpha,n,neq,q,n_post) 
    do j = 1, ny
        do i = 1, nx
            
    n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(1) = g(1,i,j)-g(3,i,j)
    n(2) = g(2,i,j)-g(4,i,j)
    n(3) = -4.0d0*g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)
        
            neq(0) = T(i,j)
            neq(1) = T(i,j)*u(i,j)
            neq(2) = T(i,j)*v(i,j)
            neq(3) = T(i,j)*paraA
            neq(4) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qnu
            q(4) = Qnu
        
            do alpha=0,4
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j) = 0.2d0*n_post(0)-0.2d0*n_post(3)
    g_post(1,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(2,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
    g_post(3,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(4,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine collisionT
!*!*End subroutine collisionT


!*!*Begin subroutine streamingT
    subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    !$omp parallel do default(none) shared(g, g_post, ex, ey, obst) private(i, j, ip, jp, alpha)
    do j = 1, ny
        do i = 1, nx
            if(obst(i,j).EQ.0) then
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                
                    g(alpha,i,j) = g_post(alpha,ip,jp)
    
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do
    
    return
    end subroutine streamingT
!*!*End subroutine streamingT


!*!*Begin subroutine bouncebackT
    subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0,y0,q,Cd1,Cd2,Cd3,Cd4


#ifdef HorizontalWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(i)
    do i = 2, nx-1
        !Bottom side
        g(2, i, 1) = g_post(4, i, 1)

        !Top side
        g(4, i, ny) = g_post(2, i, ny)
    enddo
    !$omp end parallel do
#endif

#ifdef HorizontalWallsConstT
    !$omp parallel do default(none) shared(g,g_post) private(i)
    do i = 2, nx-1
        !Bottom side
        g(2, i, 1) = -g_post(4, i, 1)+(4.0d0+paraA)/10.0d0*Thot

        !Top side
        g(4, i, ny) = -g_post(2, i, ny)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsConstT
    !$omp parallel do default(none) shared(g,g_post) private(j)
    do j = 2, ny-1
        !Left side
        g(1,1,j) = -g_post(3,1,j)+(4.0d0+paraA)/10.0d0*Thot

        !Right side
        g(3,nx,j) = -g_post(1,nx,j)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post) private(j)
    do j = 2, ny-1
        !Left side
        g(1,1,j) = g_post(3,1,j)
        
        !Right side
        g(3,nx,j) = g_post(1,nx,j)
    enddo
    !$omp end parallel do
#endif

#ifndef VerticalWallsPeriodicalT
    i = 1
    j = 1
    g(1, i, j) = -g_post(3, i, j)+(4.0d0+paraA)/10.0d0*Thot
    g(2, i, j) = -g_post(4, i, j)+(4.0d0+paraA)/10.0d0*Thot

    i = nx
    j = 1
    g(3, i, j) = -g_post(1, i, j)+(4.0d0+paraA)/10.0d0*Thot
    g(2, i, j) = -g_post(4, i, j)+(4.0d0+paraA)/10.0d0*Thot
    
    i = 1
    j = ny
    g(1, i, j) = -g_post(3, i, j)+(4.0d0+paraA)/10.0d0*Tcold
    g(4, i, j) = -g_post(2, i, j)+(4.0d0+paraA)/10.0d0*Tcold
    
    i = nx
    j = ny
    g(3, i, j) = -g_post(1, i, j)+(4.0d0+paraA)/10.0d0*Tcold
    g(4, i, j) = -g_post(2, i, j)+(4.0d0+paraA)/10.0d0*Tcold
#endif
    
#ifdef VerticalWallsPeriodicalT
    !$omp parallel do default(none) shared(g,g_post) private(j) 
    do j=1,ny
        !Left side
        g(1,1,j) = g_post(1,nx,j)

        !Right side
        g(3,nx,j) = g_post(3,1,j)
    enddo
    !$omp end parallel do
    i = 1
    j = 1
    g(2, i, j) = g_post(2, nx, j)
    
    i = nx
    j = 1
    g(2, i, j) = g_post(2, 1, j)
    
    i = 1
    j = ny
    g(4, i, j) = g_post(4, nx, j)
    
    i = nx
    j = ny
    g(4, i, j) = g_post(4, 1, j)
#endif

    return
    end subroutine bouncebackT
!*!*end subroutine bouncebackT


!*!*Begin subroutine macroT
    subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: a
    
    !$omp parallel do default(none) shared(g, T, obst, a) private(i,j)
    do j = 1, ny
        do i = 1, nx
            if(obst(i,j).EQ.0) then
                T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            endif
        enddo
    enddo
    !$omp end parallel do

    return
    end subroutine macroT
!*!*End subroutine macroT


!*!*Begin subroutine check
#ifdef steadyFlow
    subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,T,Tp,obst) private(i,j) reduction(+:error1,error2,error5,error6)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                error5 = error5+dABS( T(i,j)-Tp(i,j) )
                error6 = error6+dABS( T(i,j) )
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j)
                Tp(i,j) = T(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do 
    
    errorU = dsqrt(error1)/dsqrt(error2)
    errorT = error5/error6

    open(unit=01,file='convergence.log',status='unknown',position='append')
    write(01,*) itc,' ',errorU,' ',errorT
    close(01)
    
    return
    end subroutine check
#endif
!*!*end subroutine check

    
!*!*Begin subroutine output_binary
    subroutine output_binary()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc
#endif
#ifdef unsteadyFlow
    binFileNum = binFileNum+1
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif
    filename = adjustl(filename)

    open(unit=03,file=trim(binFolderPrefix)//"-"//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(03) ((u(i,j),i=1,nx),j=1,ny)
    write(03) ((v(i,j),i=1,nx),j=1,ny)
    write(03) ((T(i,j),i=1,nx),j=1,ny)
    close(03)

    return
    end subroutine output_binary
!*!*End subroutine output_binary

    
!*!*Begin subroutine backupData
    subroutine backupData()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,*) itc
#endif
#ifdef unsteadyFlow
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif
    filename = adjustl(filename)

    open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(05) (((f(alpha,i,j),alpha=0,8), i=1,nx), j=1,ny)
    write(05) (((g(alpha,i,j),alpha=0,4), i=1,nx), j=1,ny)
    write(05) ((u(i,j), i=1,nx), j=1,ny)
    write(05) ((v(i,j), i=1,nx), j=1,ny)
    write(05) ((T(i,j), i=1,nx), j=1,ny)
    close(05)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Backup  f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
    close(00)

    return
    end subroutine backupData
!*!*End subroutine backupData


!*!*Begin subroutine output_Tecplot
    subroutine output_Tecplot()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    character(len=100) :: filename

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1
    write(filename,'(i12.12)') pltFileNum
#endif
    filename = adjustl(filename)
    
    open(unit=41,file=trim(pltFolderPrefix)//"-"//trim(filename)//'.plt',form='binary')

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

    !c-- Number of variables in this data file

    write(41) 6

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)
    V6='rho'
    call dumpstring(V6)


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
                write(41) real(T(i,j))
                write(41) real(rho(i,j))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
    end subroutine output_Tecplot
!*!*end subroutine output_Tecplot
!*!*Begin subroutine dumpstring
    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring
!*!*end subroutine dumpstring



!*!*Begin subroutine calNuRe
    subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: halfTime
    real(kind=8) :: NuVolAvg_firstHalf, NuVolAvg_lastHalf
    real(kind=8) :: ReVolAvg_firstHalf, ReVolAvg_lastHalf
    real(kind=8) :: NuVolAvg_temp
    real(kind=8) :: ReVolAvg_temp
    real(kind=8) :: angularMomentum

    dimensionlessTime = dimensionlessTime+1
    !-------------------------------------------------------------------------------------------------------
    angularMomentum = 0.0d0
    !$omp parallel do default(none) shared(u,v,obst) private(i,j) reduction(+:angularMomentum)
    do j=1,ny 
        do i=1,nx 
            if(obst(i,j).EQ.0) then
                angularMomentum = angularMomentum+(i-nxHalf)*v(i,j)-(j-nyHalf)*u(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do
    angularMomentum = angularMomentum/dble(fluidNumMax)
    
    open(unit=01,file="angularMomentum.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), angularMomentum
    close(01)
    
    NuVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(v,T,obst) private(i,j) reduction(+:NuVolAvg_temp)
    do j=1,ny
        do i=1,nx       
            if(obst(i,j).EQ.0) then
                NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(fluidNumMax)*lengthUnit/diffusivity+1.0d0
    
    open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), NuVolAvg(dimensionlessTime)
    close(01)
    
    NuVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)+NuVolAvg(i) 
    enddo
    NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------
    ReVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(u,v,obst) private(i,j) reduction(+:ReVolAvg_temp)
    do j=1,ny
        do i=1,nx       
            if(obst(i,j).EQ.0) then
                ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            endif
        enddo
    enddo
    !$omp end parallel do
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(fluidNumMax))*lengthUnit/viscosity
    
    open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
    write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)

    ReVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)+ReVolAvg(i) 
    enddo
    ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------    

    if(statisticallyStationaryState.EQ.1) then
        if( ((dimensionlessTime-statisticallyStationaryTime).GE.minAvgPeriod).AND.(MOD(dimensionlessTime-statisticallyStationaryTime,stationaryCheckInterval).EQ.0) ) then
            
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
            
            if( (dabs((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf).LE.statConvError).AND. &
                 (dabs((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf).LE.statConvError) ) then

                if(flowReversal.EQ.0) then
                    statisticallyConvergeState = 1
                elseif(flowReversal.EQ.1) then
                    if(dimensionlessTime.GE.(statisticallyStationaryTime+flowReversalTimeMax)) statisticallyConvergeState = 1
                endif
                
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|---------Statistically Converge Reached at Time =", int(outputFrequency*dimensionlessTime),"----------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif

            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
            write(01,*) "statisticallyStationaryTime=", int(statisticallyStationaryTime*outputFrequency)
            write(01,*) "halfTime =", int(halfTime*outputFrequency)
            write(01,*) "NuVolAvg diff. (first and last half):",dabs(((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (first and last half):",dabs(((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "               "
            close(01)        
        endif
    endif

    if(statisticallyStationaryState.EQ.0) then
        if(MOD(dimensionlessTime,stationaryCheckInterval).EQ.0) then
            if( (dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/NuVolAvg_mean(dimensionlessTime)).LE.statStationError).AND. &
                (dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/ReVolAvg_mean(dimensionlessTime)).LE.statStationError) ) then
                statisticallyStationaryState = 1
                statisticallyStationaryTime = dimensionlessTime
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|-----------Statistically Stationary Reached at Time =", int(statisticallyStationaryTime*outputFrequency),"------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "calling backupData..."
                call backupData()
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif 
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
            write(01,*) "NuVolAvg diff. (stationaryCheckInterval):",dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
                                                                            /NuVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (stationaryCheckInterval):",dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
                                                                            /ReVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "               "
            close(01)       
        endif
    endif
    
    return
    end subroutine calNuRe
!*!*end subroutine calNuRe

