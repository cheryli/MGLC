!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model
!!!    Copyright (C) 2013 - 2021 Ao Xu
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

!~!~For steadyFlow, we define array up, vp, Tp to check  convergence;
!~!~for unsteadyFlow, we did not define array up, vp, Tp to save memory.
#define steadyFlow    
!~#define unsteadyFlow

!~~velocity B.C.~~
#define HorizontalWallsNoslip
!~#define VerticalWallsNoslip
#define VerticalWallsPeriodicalU
!~~velocity B.C.~~

!~~temperature B.C. (for Rayleigh Benard Cell)~~
#define RayleighBenardCell
#define HorizontalWallsConstT
!~#define VerticalWallsAdiabatic
#define VerticalWallsPeriodicalT
!~~temperature B.C.~~

!~~temperature B.C. (for Side Heated Cell)~~
!~#define SideHeatedCell
!~#define HorizontalWallsAdiabatic
!~#define VerticalWallsConstT
!~~temperature B.C.~~

!~!~For convection cells other than a simple rectangular geometry
!~#define irregularGeo
!~!~For porous convection cells, in which the solid and fluid have same thermal diffusivity
!~#define equalThermalProperty

    module ioFolder
        character(len=100) :: binFolderPrefix="../binFile/buoyancyCavity"
        character(len=100) :: pltFolderPrefix="../pltFile/buoyancyCavity"
        character(len=100) :: porousGeoFile="../case1/obstSquareCylinder-126.bin"
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
        
        integer(kind=4), parameter :: nx=513, ny=257     !----Section 1----!
        integer(kind=4), parameter :: flowGeometry=1    !----Section 1----!
        real(kind=8), parameter :: lengthUnit=dble(nx)     !----Section 1----!
        
        real(kind=8), parameter :: Rayleigh=1e5        !----Section 2----!
        real(kind=8), parameter :: Prandtl=0.71d0       !----Section 2----!
        real(kind=8), parameter :: Mach=0.1d0           !----Section 2----!
        
        real(kind=8), parameter :: outputFrequency=1.0d0 !~unit free fall time                            !----Section 3----!
        
        integer(kind=4), parameter :: dimensionlessTimeMax=int(12000/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: flowReversalTimeMax=int(10000/outputFrequency)    !----Section 3----!
                                                        !~!~if flowReversal=1,dimensionlessTimeMax shoud > flowReversalTimeMax
        integer(kind=4), parameter :: backupInterval=1000        !~!~unit: free-fall time            !----Section 3----!
        integer(kind=4), parameter :: minAvgPeriod=int(1000/outputFrequency)               !----Section 3----!
        integer(kind=4), parameter :: stationaryCheckInterval=int(200/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: convergeCheckInterval=int(200/outputFrequency)  !----Section 3----!
        
        real(kind=8), parameter :: statStationError=0.01d0  !----Section 3----!
        real(kind=8), parameter :: statConvError=0.01d0     !----Section 3----!
        real(kind=8), parameter :: epsU=1e-9                     !----Section 3----!
        real(kind=8), parameter :: epsT=1e-9                     !----Section 3----!
        
        integer(kind=4), parameter :: outputBinFile=1, outputPltFile=0                              !----Section 4----!
        !-----------------------------------------------------------------------------------------------
        !----Section 1----!
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        integer(kind=4), parameter :: nxFourth=(nx-1)/4+1, nyFourth=(ny-1)/4+1
        integer(kind=4), parameter :: nxThreeFourths=3*(nx-1)/4+1, nyThreeFourths=3*(ny-1)/4+1
        real(kind=8), parameter :: xCenter=dble(nxHalf), yCenter=dble(nyHalf)
        
        integer(kind=4), parameter :: squareGeo=1, cornerLessGeo=2, porousIsoGeo=3
        integer(kind=4), parameter :: gapX=2, gapY=3, thickness=4
        !-----------------------------------------------------------------------------------------------
        !----Section 2----!
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
        !----Section 5----!
        integer(kind=4) :: itc
        integer(kind=4), parameter :: itc_max=dimensionlessTimeMax*int(outputFrequency*timeUnit)
        
        integer(kind=4) :: binFileNum, pltFileNum
        integer(kind=4) :: dimensionlessTime
        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)  !*here, the index should start from zero
        real(kind=8) :: NuVolAvg_mean(0:dimensionlessTimeMax), ReVolAvg_mean(0:dimensionlessTimeMax)  !*because we will calculate Nu(t)-Nu(t-100) at time t=100.
        !-----------------------------------------------------------------------------------------------
        !----Section 6----!
        integer(kind=4) :: fluidNumMax
        integer(kind=4), allocatable :: obst(:,:)
        !-----------------------------------------------------------------------------------------------     
    end module commondata


!*!*Begin program main
    program main
    use commondata
    use ioFolder
    implicit none
    real(kind=8) :: timeStart, timeEnd
    character(len=24) :: ctime, string
    INTEGER(kind=4) :: time
    
    open(unit=00,file="SimulationSettings.txt",status='unknown')
    string = ctime( time() )
    write(00,*) 'Start: ', string
    close(00)

    call initial()

    call CPU_TIME(timeStart)

#ifdef steadyFlow
    !$acc data copy(u,v,T) copyin(f,g,up,vp,Tp,rho,obst,omega,omegaT,ex,ey,r) create(f_post,g_post,Fx,Fy)
#endif
#ifdef unsteadyFlow
    !$acc data copy(u,v,T) copyin(f,g,rho,obst,omega,omegaT,ex,ey,r) create(f_post,g_post,Fx,Fy)
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
        
        if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then
            
#ifdef steadyFlow
            if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()
#endif
            
#ifdef unsteadyFlow
            if(statisticallyStationaryState.EQ.1) then
                if(outputBinFile.EQ.1) then
                    call output_binary()
                    if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()
                endif   
                if(outputPltFile.EQ.1) then
                    call output_Tecplot()
                endif
            endif
#endif
            
            call calNuRe()
            
        endif

    enddo
    
    if(outputBinFile.EQ.1) call backupData()
    
    if(outputPltFile.EQ.1) call output_Tecplot()
    !$acc end data

    call CPU_TIME(timeEnd)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Time (CPU) = ", real(timeEnd-timeStart), "s"

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
        write(00,*) "----------------------------------"
        write(00,*) "paraA=", paraA
        write(00,*) "Error: condition not meet for the algorithm"
        write(00,*) "Ref: Luo2013, CMA"
        write(00,*) "Please try to reduce Mach number"
        write(00,*) "----------------------------------"
        stop
    endif
    
    if(flowReversal.EQ.1) then
        if(dimensionlessTimeMax.LE.flowReversalTimeMax) then
            write(00,*) "Error: in case of flor reversal,"
            write(00,*) "...please check dimensionlessTimeMax!"
            stop
        endif
    endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',nx,ny
    write(00,*) 'Rayleigh=',real(Rayleigh), '; Prandtl =',real(Prandtl), '; Mach =',real(Mach)
    write(00,*) "Length unit: L0 =", real(lengthUnit)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf)
    write(00,*) "viscosity =",real(viscosity), "; diffusivity =",real(diffusivity)
    write(00,*) "outputFrequency =", real(outputFrequency), "tf"
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

#ifdef irregularGeo
    write(00,*) "I am irregular geometry"
    if(flowGeometry.EQ.cornerLessGeo) then
        write(00,*) "I am cornerLessGeo"
        write(00,*) "gapX=",gapX,"; gapY=",gapY,"; thickness=",thickness
    elseif(flowGeometry.EQ.porousIsoGeo) then
        write(00,*) "I am porousIsoGeo"
        write(00,*) "porous geometry information load from external file", porousGeoFile
    elseif(flowGeometry.EQ.squareGeo) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        stop
    endif
#endif
#ifndef irregularGeo
    write(00,*) "I am regular geometry (Square Cavity)"
    if(flowGeometry.NE.squareGeo) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "flowGeometry should be", squareGeo, "; while its actual value is ", flowGeometry
        stop
    endif
    if( (dabs(lengthUnit-dble(nx)).GT.1e-6).AND.(dabs(lengthUnit-dble(ny)).GT.1e-6) ) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "Tip: lengthUnit should be equal to nx or ny"
        stop
    endif
#endif

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

    allocate (f(nx,ny,0:8))
    allocate (f_post(0:nx+1,0:ny+1,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(0:nx+1,0:ny+1,0:4))

    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    allocate (obst(0:nx+1,0:ny+1))
    
    obst = 0
    rho = 1.0d0
    
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
    
#ifdef irregularGeo
    if(flowGeometry.EQ.cornerLessGeo) then
        fluidNumMax = nx*ny
        do j=1+gapY,ny-gapY
            do i=1+gapX,nxFourth+gapX
                if( ((i+j).GE.(nyFourth+gapX+gapY)).AND.((i+j).LE.(nyFourth+gapX+gapY+thickness)) ) then
                    obst(i,j) = 1 !~!~Bottom left corner
                    fluidNumMax = fluidNumMax-1
                endif
                if( ((-i+j).GE.(nyThreeFourths-gapX-gapY-thickness)).AND.((-i+j).LE.(nyThreeFourths-gapX-gapY)) ) then
                    obst(i,j) = 1 !~!~ Top left corner
                    fluidNumMax = fluidNumMax-1
                endif
            enddo
        enddo
        do j=1+gapY,ny-gapY
            do i=nx-nxFourth+1-gapX,nx-gapX
                if( ((-i+j).GE.(-nyThreeFourths+gapX+gapY)).AND.((-i+j).LE.(-nyThreeFourths+gapX+gapY+thickness)) ) then
                    obst(i,j) = 1 !~!~Bottom right corner
                    fluidNumMax = fluidNumMax-1
                endif
                if( ((i+j).GE.(ny+nyThreeFourths-gapX-gapY-thickness+1)).AND.((i+j).LE.(ny+nyThreeFourths-gapX-gapY+1)) ) then
                    obst(i,j) = 1 !~!~Top right corner
                    fluidNumMax = fluidNumMax-1
                endif
            enddo
        enddo
    elseif(flowGeometry.EQ.porousIsoGeo) then
        open(unit=01,file=porousGeoFile,form="unformatted",access="sequential",status='old')
        read(01) ((obst(i,j),i=1,nx),j=1,ny)
        close(01)
        fluidNumMax = nx*ny
        do j=1,ny
            do i=1,nx
                if(obst(i,j).EQ.1) fluidNumMax=fluidNumMax-1
            enddo
        enddo
    endif
#endif

    if(loadInitField.EQ.0) then 
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
        
        write(00,*) "Initial field is set exactly"
        if(reloadDimensionlessTime.NE.0) then
            write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
            stop
        endif
        
#ifdef VerticalWallsNoslip
        write(00,*) "Velocity B.C. for vertical walls are: ===No-slip wall==="
#endif
#ifdef VerticalWallsPeriodicalU
        write(00,*) "Velocity B.C. for vertical walls are: ===Periodical==="
#endif
#ifdef HorizontalWallsNoslip
        write(00,*) "Velocity B.C. for horizontal walls are: ===No-slip wall==="
#endif

#ifdef irregularGeo
        write(00,*) "Velocity B.C. for irregular curved walls are: ===No-slip wall==="
        if(flowGeometry.EQ.cornerLessGeo) write(00,*)"Temperature B.C. for irregular curved walls are:===Adiabatic wall==="
        if(flowGeometry.EQ.porousIsoGeo) write(00,*) "Temperature B.C. for irregular curved walls are:===Adiabatic wall==="
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

#ifdef VerticalWallsAdiabatic
    write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
    write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif

        f = 0.0d0
        g = 0.0d0
        do j=1,ny
            do i=1,nx
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha=0,8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(i,j,alpha) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo
                do alpha=0,4
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    g(i,j,alpha) = omegaT(alpha)*T(i,j)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))
                enddo
            enddo
        enddo

    elseif(loadInitField.EQ.1) then
        if(reloadDimensionlessTime.EQ.0) then
            write(00,*) "Error: since loadInitField.EQ.1, reloadDimensionlessTime should not be 0"
            stop
        endif
        write(00,*) "Load initial field from previous simulation: ../reloadFile/backupFile- >>>"
        write(reloadFileName, *) reloadbinFileNum
        reloadFileName = adjustl(reloadFileName)
        open(unit=01,file="../reloadFile/backupFile-"//trim(reloadFileName)//".bin",form="unformatted",access="sequential",status='old')
        read(01) (((f(i,j,alpha), i=1,nx), j=1,ny), alpha=0,8)
        read(01) (((g(i,j,alpha), i=1,nx), j=1,ny), alpha=0,4)
        read(01) ((u(i,j), i=1,nx), j=1,ny)
        read(01) ((v(i,j), i=1,nx), j=1,ny)
        read(01) ((T(i,j), i=1,nx), j=1,ny)
        close(01)

        dev = 0.0d0
        do j=1,ny
            do i=1,nx
                rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
                dev = dev+g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)-T(i,j)
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
    dimensionlessTime = 0
    
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 0
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

    !$acc parallel loop &
    !$acc private(m,meq,m_post,s,fSource) present(f,f_post,rho,u,v,T,Fx,Fy,obst,omega) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
        
#ifdef irregularGeo
            if(obst(i,j).EQ.0) then
#endif

    m(0) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
    m(1) = -4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+2.0d0*(f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8))
    m(2) = 4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
    m(3) = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
    m(4) = -2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
    m(5) = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
    m(6) = -2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
    m(7) = f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
    m(8) = f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

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

    f_post(i,j,0) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0
    f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
    f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

#ifdef irregularGeo
            else
                do alpha=0,8
                    f(i,j,alpha) = rho(i,j)*omega(alpha)
                    f_post(i,j,alpha) = rho(i,j)*omega(alpha)
                enddo
            endif
#endif

        enddo
    enddo
    !$acc end parallel loop
    
    return
    end subroutine collision
!*!*End subroutine collision


!*!*Begin subroutine streaming
    subroutine streaming()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    
    !$acc parallel loop present(f,f_post,obst) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
        
            if(obst(i,j).EQ.0) then
            
                f(i,j,0) = f_post(i,j,0)
                
                f(i,j,1) = f_post(i-1,j,1)
                f(i,j,2) = f_post(i,j-1,2)
                f(i,j,3) = f_post(i+1,j,3)
                f(i,j,4) = f_post(i,j+1,4)
                
                f(i,j,5) = f_post(i-1,j-1,5)
                f(i,j,6) = f_post(i+1,j-1,6)
                f(i,j,7) = f_post(i+1,j+1,7)
                f(i,j,8) = f_post(i-1,j+1,8)
            
            endif
                
        enddo
    enddo
    !$acc end parallel loop
    
    return
    end subroutine streaming
!*!*End subroutine streaming


!*!*Begin subroutine bounceback
    subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    
#ifdef irregularGeo
    !$acc parallel loop present(f,f_post,obst,ex,ey,r) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=1,8
                    if(obst(i+ex(alpha),j+ey(alpha)).EQ.1) then
                        if(flowGeometry.EQ.porousIsoGeo) then
                            f(i,j,r(alpha)) = f_post(i,j,alpha)
                        endif
                    endif
                enddo
            endif
        enddo
    enddo
#endif

#ifdef VerticalWallsPeriodicalU
    !$acc parallel loop present(f,f_post)
    do j=1,ny
        !Left side (i=1)
        f(1,j,1) = f_post(nx,j,1)
        f(1,j,5) = f_post(nx,j,5)
        f(1,j,8) = f_post(nx,j,8)
        
        !Right side (i=nx)
        f(nx,j,3) = f_post(1,j,3)
        f(nx,j,6) = f_post(1,j,6)
        f(nx,j,7) = f_post(1,j,7)
    enddo
    !$acc end parallel loop
#endif

#ifdef VerticalWallsNoslip
    !$acc parallel loop present(f,f_post)
    do j=1,ny
        !Left side (i=1)
        f(1,j,1) = f_post(1,j,3)
        f(1,j,5) = f_post(1,j,7)
        f(1,j,8) = f_post(1,j,6)

        !Right side (i=nx)
        f(nx,j,3) = f_post(nx,j,1)
        f(nx,j,6) = f_post(nx,j,8)
        f(nx,j,7) = f_post(nx,j,5)
    enddo
    !$acc end parallel loop
#endif

#ifdef HorizontalWallsNoslip
    !$acc parallel loop present(f,f_post) 
    do i=1,nx
        !Bottom side (j=1)
        f(i,1,2) = f_post(i,1,4)
        f(i,1,5) = f_post(i,1,7)
        f(i,1,6) = f_post(i,1,8)

        !Top side (j=ny)
        f(i,ny,4) = f_post(i,ny,2)
        f(i,ny,7) = f_post(i,ny,5)
        f(i,ny,8) = f_post(i,ny,6)
    enddo
    !$acc end parallel loop
#endif

    return
    end subroutine bounceback
!*!*End subroutine bounceback


!*!*Begin subroutine macro
    subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop present(f,rho,u,v,Fx,Fy,obst) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
            if(obst(i,j).EQ.0) then
                rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
                u(i,j) = ( f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*Fx(i,j) )/rho(i,j)
                v(i,j) = ( f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*Fy(i,j) )/rho(i,j)
            elseif(obst(i,j).EQ.1) then
                rho(i,j) = 1.0d0
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0
            endif
        enddo
    enddo
    !$acc end parallel loop
    
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

    !$acc parallel loop &
    !$acc private(n,neq,q,n_post) present(g,g_post,u,v,T,obst,omegaT)
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
            
#ifdef irregularGeo
#ifndef equalThermalProperty
            if(obst(i,j).EQ.0) then
#endif
#endif

    n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
    n(1) = g(i,j,1)-g(i,j,3)
    n(2) = g(i,j,2)-g(i,j,4)
    n(3) = -4.0d0*g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
    n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)
        
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
        
    g_post(i,j,0) = 0.2d0*n_post(0)-0.2d0*n_post(3)
    g_post(i,j,1) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(i,j,2) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
    g_post(i,j,3) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(i,j,4) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        
#ifdef irregularGeo
#ifndef equalThermalProperty
            else
                do alpha=0,4
                    g(i,j,alpha) = T(i,j)*omegaT(alpha)
                    g_post(i,j,alpha) = T(i,j)*omegaT(alpha)
                enddo
            endif
#endif
#endif

        enddo
    enddo
    !$acc end parallel loop
    
    return
    end subroutine collisionT
!*!*End subroutine collisionT


!*!*Begin subroutine streamingT
    subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop present(g,g_post,obst) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
        
#ifndef equalThermalProperty
            if(obst(i,j).EQ.0) then
#endif
                g(i,j,0) = g_post(i,j,0)

                g(i,j,1) = g_post(i-1,j,1)
                g(i,j,2) = g_post(i,j-1,2)
                g(i,j,3) = g_post(i+1,j,3)
                g(i,j,4) = g_post(i,j+1,4)
            
#ifndef equalThermalProperty
            endif
#endif
            
        enddo
    enddo
    !$acc end parallel loop

    return
    end subroutine streamingT
!*!*End subroutine streamingT


!*!*Begin subroutine bouncebackT
    subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    
#ifdef irregularGeo
#ifndef equalThermalProperty
    !$acc parallel loop present(g,g_post,obst,ex,ey,r) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=1,4
                    if(obst(i+ex(alpha),j+ey(alpha)).EQ.1) then
                        if(flowGeometry.EQ.porousIsoGeo) then
                            g(i,j,r(alpha)) = g_post(i,j,alpha)
                        endif
                    endif
                enddo
            endif
        enddo
    enddo
#endif
#endif

#ifdef HorizontalWallsAdiabatic
    !$acc parallel loop present(g,g_post)
    do i=1,nx
        !Bottom side
        g(i,1,2) = g_post(i,1,4)

        !Top side
        g(i,ny,4) = g_post(i,ny,2)
    enddo
    !$acc end parallel loop
#endif

#ifdef HorizontalWallsConstT
    !$acc parallel loop present(g,g_post)
    do i=1,nx
        !Bottom side
        g(i,1,2) = -g_post(i,1,4)+(4.0d0+paraA)/10.0d0*Thot

        !Top side
        g(i,ny,4) = -g_post(i,ny,2)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
    !$acc end parallel loop
#endif

#ifdef VerticalWallsConstT
    !$acc parallel loop present(g,g_post) 
    do j=1,ny
        !Left side
        g(1,j,1) = -g_post(1,j,3)+(4.0d0+paraA)/10.0d0*Thot

        !Right side
        g(nx,j,3) = -g_post(nx,j,1)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
    !$acc end parallel loop
#endif

#ifdef VerticalWallsAdiabatic
    !$acc parallel loop present(g,g_post)
    do j=1,ny
        !Left side
        g(1,j,1) = g_post(1,j,3)

        !Right side
        g(nx,j,3) = g_post(nx,j,1)
    enddo
    !$acc end parallel loop
#endif

#ifdef VerticalWallsPeriodicalT
    !$acc parallel loop present(g,g_post)
    do j=1,ny
        !Left side
        g(1,j,1) = g_post(nx,j,1)

        !Right side
        g(nx,j,3) = g_post(1,j,3)
    enddo
    !$acc end parallel loop
#endif

    return
    end subroutine bouncebackT
!*!*end subroutine bouncebackT


!*!*Begin subroutine macroT
    subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$acc parallel loop present(g,T,obst) 
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
#ifndef equalThermalProperty
            if(obst(i,j).EQ.0) then
#endif
                T(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
#ifndef equalThermalProperty
            endif
#endif
        enddo
    enddo
    !$acc end parallel loop
    
    return
    end subroutine macroT
!*!*end subroutine macroT


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
    
    !$acc parallel loop reduction(+:error1,error2,error5,error6) present(u,v,up,vp,T,Tp,obst)
    do j=1,ny
        !$acc loop gang vector
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
    !$acc end parallel loop
    
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
    
    !$acc update self(u, v, T, rho)
    
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
    write(03) ((rho(i,j), i=1,nx), j=1,ny)
    close(03)

    return
    end subroutine output_binary
!*!*end subroutine output_binary

    
!*!*Begin subroutine backupData
    subroutine backupData()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j, alpha
    character(len=100) :: filename
    
    !$acc update self(u, v, T, f, g)

#ifdef steadyFlow
    write(filename,*) itc
#endif
#ifdef unsteadyFlow
    if(loadInitField.EQ.0) write(filename,*) binFileNum
    if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
#endif
    filename = adjustl(filename)

    open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(05) (((f(i,j,alpha), i=1,nx), j=1,ny), alpha=0,8)
    write(05) (((g(i,j,alpha), i=1,nx), j=1,ny), alpha=0,4)
    write(05) ((u(i,j), i=1,nx), j=1,ny)
    write(05) ((v(i,j), i=1,nx), j=1,ny)
    write(05) ((T(i,j), i=1,nx), j=1,ny)
    close(05)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Backup  f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
    close(00)
    
    return
    end subroutine backupData
!*!*end subroutine backupData
    
    
!*!*Begin subroutine output_Tecplot
    subroutine output_Tecplot()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    character(len=100) :: filename
    
    !$acc update self(u,v,T)

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1
    write(filename,'(i12.12)') pltFileNum
#endif
    filename = adjustl(filename)

    open(41,file=trim(pltFolderPrefix)//trim(filename)//'.plt', access='stream', form='unformatted')

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

    !c-- Number of variables in this data file
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
    V5='T'
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
                write(41) real(T(i,j))
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
    !$acc parallel loop reduction(+:angularMomentum) present(u,v,obst)
    do j=1,ny 
        !$acc loop gang vector
        do i=1,nx 
            if(obst(i,j).EQ.0) then
                angularMomentum = angularMomentum+(i-nxHalf)*v(i,j)-(j-nyHalf)*u(i,j)
            endif
        enddo
    enddo
    !$acc end parallel loop
    angularMomentum = angularMomentum/dble(fluidNumMax)
    
    open(unit=01,file="angularMomentum.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), angularMomentum
    close(01)
    
    NuVolAvg_temp = 0.0d0    
    !$acc parallel loop reduction(+:NuVolAvg_temp) present(v,T,obst)
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx
            if(obst(i,j).EQ.0) then
                NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)
            endif
        enddo
    enddo
    !$acc end parallel loop
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
    !$acc parallel loop reduction(+:ReVolAvg_temp) present(u,v,obst)
    do j=1,ny
        !$acc loop gang vector
        do i=1,nx 
            if(obst(i,j).EQ.0) then
                ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            endif
        enddo
    enddo
    !$acc end parallel loop
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
        if( ((dimensionlessTime-statisticallyStationaryTime).GE.minAvgPeriod).AND.(MOD(dimensionlessTime-statisticallyStationaryTime,convergeCheckInterval).EQ.0) ) then
            
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
                write(01,*) "|---------Statistically Converge Reached at Time =", int(dimensionlessTime*outputFrequency),"----------|"
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
        if(MOD(dimensionlessTime, stationaryCheckInterval).EQ.0) then
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

