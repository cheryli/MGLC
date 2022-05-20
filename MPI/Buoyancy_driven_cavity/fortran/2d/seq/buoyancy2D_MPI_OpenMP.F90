!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model
!!!    Copyright (C) 2013 - 2021  Ao Xu
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

!~!~For steadyFlow, we define array up, vp, Tp to check  convergence;
!~!~for unsteadyFlow, we did not define array up, vp, Tp to save memory.
#define steadyFlow    
!~#define unsteadyFlow

!~~velocity B.C.~~
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!~#define VerticalWallsPeriodicalU
!~~velocity B.C.~~

!~~temperature B.C. (for Rayleigh Benard Cell)~~
! #define RayleighBenardCell
!~#define HorizontalWallsConstT
!~#define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT
!~~temperature B.C.~~

!~~temperature B.C. (for Side Heated Cell)~~
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!~~temperature B.C.~~

!~!~For convection cells other than a simple rectangular geometry
! #define irregularGeo
!~!~For porous convection cells, in which the solid and fluid have same thermal diffusivity
!~#define equalThermalProperty

    module ioFolder
        character(len=30) :: binFolderPrefix="../binFile/buoyancyCavity"
        character(len=100) :: pltFolderPrefix="../pltFile/buoyancyCavity"
        character(len=100) :: porousGeoFile="../case1/obstSquareCylinder-126.bin"
    end module ioFolder

    module commondata
        use mpi 
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
        
        integer(kind=4), parameter :: nx=2049, ny=2049     !----Section 1----!
        integer(kind=4), parameter :: flowGeometry=1    !----Section 1----!
        integer(kind=4), parameter :: heaterLength=1 !! if flowGeometry.EQ.4, then this value is in use
        real(kind=8), parameter :: margin=8.0d0 !! should be an even number, may not  in use  
        real(kind=8), parameter :: diameter0=dble(ny)-margin !! if flowGeometry.EQ.4, then this value is in use
        ! if flowGeometry.EQ.1, lengthUnit.EQ.ny
        ! if flowGeometry.EQ.4, lengthUnit.EQ.diameter0
        real(kind=8), parameter :: lengthUnit=dble(ny)     !----Section 1----!  
        
        real(kind=8), parameter :: Rayleigh=1e7        !----Section 2----!
        real(kind=8), parameter :: Prandtl=0.71d0       !----Section 2----!
        real(kind=8), parameter :: Mach=0.1d0           !----Section 2----!
        
        real(kind=8), parameter :: outputFrequency=1.0d0 !~unit free fall time                            !----Section 3----!
        integer(kind=4), parameter :: dimensionlessTimeMax=int(5000/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: flowReversalTimeMax=int(200/outputFrequency)    !----Section 3----!
                                                        !~!~if flowReversal=1,dimensionlessTimeMax shoud > flowReversalTimeMax
        integer(kind=4), parameter :: backupInterval=2000        !~!~unit: free-fall time            !----Section 3----!
        integer(kind=4), parameter :: minAvgPeriod=int(1000/outputFrequency)               !----Section 3----!
        integer(kind=4), parameter :: stationaryCheckInterval=int(200/outputFrequency)  !----Section 3----!
        integer(kind=4), parameter :: convergeCheckInterval=int(200/outputFrequency)  !----Section 3----!
        
        real(kind=8), parameter :: statStationError=0.01d0  !----Section 3----!
        real(kind=8), parameter :: statConvError=0.01d0     !----Section 3----!
        real(kind=8), parameter :: epsU=1e-6                     !----Section 3----!
        real(kind=8), parameter :: epsT=1e-6                     !----Section 3----!
        
        integer(kind=4), parameter :: outputBinFile=0, outputPltFile=0                              !----Section 4----!
        !-----------------------------------------------------------------------------------------------
        !----Section 1----!
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        integer(kind=4), parameter :: nxFourth=(nx-1)/4+1, nyFourth=(ny-1)/4+1
        integer(kind=4), parameter :: nxThreeFourths=3*(nx-1)/4+1, nyThreeFourths=3*(ny-1)/4+1
        real(kind=8), parameter :: xCenter=dble(nxHalf), yCenter=dble(nyHalf)
        integer(kind=4), parameter :: squareGeo=1, cornerLessGeo=2, porousIsoGeo=3, circularGeo=4
        integer(kind=4), parameter :: oneThirdHeater=3, oneFourthHeater=4, allHeater=1
        real(kind=8), parameter :: radius0=diameter0/2.0d0
        integer(kind=4), parameter :: gapX=2, gapY=3, thickness=4
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
        integer(kind=4) :: dimensionlessTime
        integer(kind=4), parameter :: itc_max=dimensionlessTimeMax*int(outputFrequency*timeUnit)
        real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)  !*here, the index should start from zero
        real(kind=8) :: NuVolAvg_mean(0:dimensionlessTimeMax), ReVolAvg_mean(0:dimensionlessTimeMax)  !*because we will calculate Nu(t)-Nu(t-100) at time t=100.

        integer(kind=4) :: binFileNum, pltFileNum
        !-----------------------------------------------------------------------------------------------
        !----Section 6----!
        integer(kind=4) :: fluidNumMax
        integer(kind=4), allocatable :: obst(:,:), obstFit(:,:)
        !-----------------------------------------------------------------------------------------------  
        !----Section 7 (for MPI)----!        
        integer(kind=4) :: nProc
        integer(kind=4), allocatable :: start1d(:), end1d(:), count1d(:), displ1d(:)
        integer(kind=4), allocatable :: start2d(:), end2d(:), count2d(:), displ2d(:)
        integer(kind=4), allocatable :: startD2Q9(:), endD2Q9(:), countD2Q9(:), displD2Q9(:)
        integer(kind=4), allocatable :: startD2Q5(:), endD2Q5(:), countD2Q5(:), displD2Q5(:)
        integer(kind=4) :: myID, iStatus(MPI_Status_Size)
        integer(kind=4) :: idest, ierr, iRoot
        integer(kind=4) :: iSrc, iTag
        integer(kind=4) :: iStart, iEnd
        integer(kind=4) :: iStartMinus1, iEndPlus1
        integer(kind=4) :: leftNeighbor, rightNeighbor
        integer(kind=4) :: nyLocal
        real(kind=8), allocatable :: uGlobal(:,:), vGlobal(:,:), tGlobal(:,:)
        real(kind=8), allocatable :: rhoGlobal(:,:)
        real(kind=8), allocatable :: fGlobal(:,:,:), gGlobal(:,:,:)
        integer(kind=4), allocatable :: obstGlobal(:,:)
        !-----------------------------------------------------------------------------------------------     
    end module commondata
    
    
!*!*Begin program main
    program main
    use mpi
    use omp_lib    
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    integer(kind=4) :: myMaxThreads
    character(len=24) :: ctime, string
    integer(kind=4) :: time
    
    !----------------
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)

    if (myID.EQ.0) then
        open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
        string = ctime( time() )
        write(00,*) 'Start: ', string
        write(00,*) "Max Running processes:", nProc
        close(00)
    endif
    
#ifdef _OPENMP
    call OMP_set_num_threads(12)
    myMaxThreads = OMP_get_max_threads()
    if (myID.EQ.0) then
        open(unit=00, file="SimulationSettings.txt", status="unknown", position="append")
        write(00,*) "Starting OpenMP >>>>>>"
        write(00,*) "Max Running threads:",myMaxThreads
        close(00)
    endif
#endif

    allocate (start1d(0:nProc-1))
    allocate (end1d(0:nProc-1))
    allocate (count1d(0:nProc-1))
    allocate (displ1d(0:nProc-1))
    
    allocate (start2d(0:nProc-1))
    allocate (end2d(0:nProc-1))
    allocate (count2d(0:nProc-1))
    allocate (displ2d(0:nProc-1))
    
    allocate (startD2Q9(0:nProc-1))
    allocate (endD2Q9(0:nProc-1))
    allocate (countD2Q9(0:nProc-1))
    allocate (displD2Q9(0:nProc-1))
    
    allocate (startD2Q5(0:nProc-1))
    allocate (endD2Q5(0:nProc-1))
    allocate (countD2Q5(0:nProc-1))
    allocate (displD2Q5(0:nProc-1))
    !--------------------------------------------------------------------
    
    nyLocal = ny/nProc+1
    allocate (u(nx,nyLocal))
    allocate (v(nx,nyLocal))
    allocate (rho(nx,nyLocal))
    allocate (T(nx,nyLocal))
#ifdef steadyFlow
    allocate (up(nx,nyLocal))
    allocate (vp(nx,nyLocal))
    allocate (Tp(nx,nyLocal))
#endif
    allocate (Fx(nx,nyLocal))
    allocate (Fy(nx,nyLocal))    
    allocate (f(0:8,nx,nyLocal))
    allocate (f_post(0:8,-1:nx+2,-1:nyLocal+2))
    allocate (g(0:4,nx,nyLocal))
    allocate (g_post(0:4,-1:nx+2,-1:nyLocal+2))
    allocate (obst(0:nx+1,-1:nyLocal+2))
    allocate (obstFit(nx,0:nyLocal+1))
    
    allocate (rhoGlobal(nx,ny))
    allocate (uGlobal(nx,ny))
    allocate (vGlobal(nx,ny))
    allocate (tGlobal(nx,ny))
    allocate (fGlobal(0:8,nx,ny))
    allocate (gGlobal(0:4,nx,ny))
    allocate (obstGlobal(nx,ny))
    !--------------------------------------------------------------------    
    call StartEnd(1, ny)

    iStart = 1
    iEnd = count1d(myID)
    !----------------
    
    !================================================
    iStartMinus1 = iStart-1
    iEndPlus1 = iEnd+1
    
    leftNeighbor = myID-1
    rightNeighbor = myID+1
    if(myID.EQ.0) leftNeighbor = MPI_PROC_NULL
    if(myID.EQ.nProc-1) rightNeighbor = MPI_PROC_NULL
    !================================================

    call initial()
    !----------------------------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    timeStart = MPI_WTIME()
    !----------------------------------------------------------

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

        if (itc == 100) then
            exit
        endif
        
! #ifdef RayleighBenardCell
!         if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then
!     !-------------------------------------------------------------------------------------------------------------------------------
!     iDest = 0
!     call MPI_GATHERV(u(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
!                            uGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
    
!     iDest = 0
!     call MPI_GATHERV(v(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
!                            vGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                            
!     iDest = 0
!     call MPI_GATHERV(T(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
!                            tGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                           
!     if(MOD(itc,backupInterval*int(timeUnit)).EQ.0) then
!         iDest = 0
!         call MPI_GATHERV(f(0,1,1), countD2Q9(myID), MPI_DOUBLE_PRECISION, &
!                                fGlobal, countD2Q9, displD2Q9, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                               
!         iDest = 0
!         call MPI_GATHERV(g(0,1,1), countD2Q5(myID), MPI_DOUBLE_PRECISION, &
!                                gGlobal, countD2Q5, displD2Q5, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
!     endif
!     !----------------------------------------------------------
            
! if(myID.EQ.0) then
! #ifdef steadyFlow
!                 if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()
! #endif
                                     
!                 if(statisticallyStationaryState.EQ.1) then
!                     if(outputBinFile.EQ.1) then
!                         call output_binary()
!                         if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()
!                     endif
!                     if(outputPltFile.EQ.1) then
!                         call output_Tecplot()
!                     endif
!                 endif
            
!                 call calNuRe()
! endif
        
!             iRoot = 0
!             call MPI_BCAST(dimensionlessTime, 1,  MPI_integer4, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(statisticallyStationaryState, 1,  MPI_integer4, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(statisticallyStationaryTime, 1,  MPI_integer4, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(statisticallyConvergeState, 1,  MPI_integer4, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(NuVolAvg, dimensionlessTimeMax+1,  MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(ReVolAvg, dimensionlessTimeMax+1,  MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(NuVolAvg_mean, dimensionlessTimeMax+1,  MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
!             call MPI_BCAST(ReVolAvg_mean, dimensionlessTimeMax+1,  MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
            
!         endif
! #endif

    enddo !~!~end of main loop

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    timeEnd = MPI_WTIME()

    if(myID.EQ.0) then
        write(*,*) "MPI time = ", real(timeEnd-timeStart),"seconds"
    endif
    
    !-------------------------------------------------------------------------------------------------------------------------------
    iDest = 0
    call MPI_GATHERV(u(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
                           uGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
    
    iDest = 0
    call MPI_GATHERV(v(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
                           vGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                            
    iDest = 0
    call MPI_GATHERV(T(1,1), count2d(myID), MPI_DOUBLE_PRECISION, &
                           tGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                           
    iDest = 0
    call MPI_GATHERV(f(0,1,1), countD2Q9(myID), MPI_DOUBLE_PRECISION, &
                           fGlobal, countD2Q9, displD2Q9, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                           
    iDest = 0
    call MPI_GATHERV(g(0,1,1), countD2Q5(myID), MPI_DOUBLE_PRECISION, &
                           gGlobal, countD2Q5, displD2Q5, MPI_DOUBLE_PRECISION, iDest, MPI_COMM_WORLD, ierr)
                           
if(myID.EQ.0) then
    if(outputBinFile.EQ.1) call backupData()
        
    if(outputPltFile.EQ.1) call output_Tecplot()
endif
    
    deallocate (uGlobal)
    deallocate (vGlobal)
    deallocate (tGlobal)
    deallocate (fGlobal)
    deallocate (gGlobal)

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
    !----------------------------------------------------------
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! timeEnd = MPI_WTIME()
! if(myID.EQ.0) then
!     open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
!     write(00,*) "   "
!     write(00,*) "Time (CPU) =", real(timeEnd-timeStart),"seconds"
!     close(00)
! endif
    !----------------------------------------------------------
    call MPI_finalize(ierr)
    
! if(myID.EQ.0) then
!     open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
!     write(00,*) "Successfully: DNS completed!"
!     string = ctime( time() )
!     write(00,*) 'End:   ', string
!     close(00)
! endif

    stop
    end program main
!*!*end program main


!*!*Begin subroutine StartEnd    
    subroutine StartEnd(iS1, iS2)
    use commondata
    implicit none
    integer(kind=4) :: leng, iBlock
    integer(kind=4) :: ir
    integer(kind=4) :: iS1, iS2
    integer(kind=4) :: i
    
    leng = iS2-iS1+1
    iBlock = leng/nProc
    ir= leng-iBlock*nProc
    
    do i=0,nProc-1
        
        if(i.LT.ir) then  
            start1d(i) = iS1+i*(iBlock+1) 
            count1d(i) = iBlock+1
            end1d(i) = start1d(i)+count1d(i)-1

            !-----------------------------------------------------------
            start2d(i) = iS1+i*(iBlock+1)*nx
            count2d(i) = (iBlock+1)*nx
            end2d(i) = start2d(i)+count2d(i)-1
            
            startD2Q9(i) = iS1+i*(iBlock+1)*nx*9
            countD2Q9(i) = (iBlock+1)*nx*9
            endD2Q9(i) = startD2Q9(i)+countD2Q9(i)-1
            
            startD2Q5(i) = iS1+i*(iBlock+1)*nx*5
            countD2Q5(i) = (iBlock+1)*nx*5
            endD2Q5(i) = startD2Q5(i)+countD2Q5(i)-1

        else
            start1d(i) = iS1+i*iBlock+ir   
            count1d(i) = iBlock   
            end1d(i) = start1d(i)+count1d(i)-1  
            
            !-----------------------------------------------------------
            start2d(i) = iS1+i*iBlock*nx+ir*nx  
            count2d(i) = iBlock*nx  
            end2d(i) = start2d(i)+count2d(i)-1  

            startD2Q9(i) = iS1+i*iBlock*nx*9+ir*nx*9
            countD2Q9(i) = iBlock*nx*9
            endD2Q9(i) = startD2Q9(i)+countD2Q9(i)-1
            
            startD2Q5(i) = iS1+i*iBlock*nx*5+ir*nx*5
            countD2Q5(i) = iBlock*nx*5
            endD2Q5(i) = startD2Q5(i)+countD2Q5(i)-1
        endif
        
        displ1d(i) = start1d(i)-iS1  
        displ2d(i) = start2d(i)-iS1  
        displD2Q9(i) = startD2Q9(i)-iS1
        displD2Q5(i) = startD2Q5(i)-iS1

    enddo
    
    return
    end subroutine StartEnd
!*!*End subroutine StartEnd
    

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
    integer(kind=4) :: tempJ
    real(kind=8) :: dev
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    fluidNumMax = nx*ny

if(myID.EQ.0) then
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
        write(00,*) 'Error, paraA=', paraA
        write(00,*) "Current Mach number is :", real(Mach)
        write(00,*) 'Please try to reduce Mach number!'
        write(00,*) '----------------------------------'
        call MPI_finalize(ierr)
        stop
    endif
        
    if(flowReversal.EQ.1) then
        if(dimensionlessTimeMax.LE.flowReversalTimeMax) then
            write(00,*) "Error: in case of flor reversal,"
            write(00,*) "...please check dimensionlessTimeMax!"
            call MPI_finalize(ierr)
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
        
#ifdef irregularGeo
    write(00,*) "I am irregular geometry"
    if(flowGeometry.EQ.circularGeo) then
        write(00,*) "I am circular cell"
        if(heaterLength.EQ.oneThirdHeater) then
            write(00,*) "I am 1/3 heater length"
        elseif(heaterLength.EQ.oneFourthHeater) then
            write(00,*) "I am 1/4 heater length"
        elseif(heaterLength.EQ.allHeater) then
            write(00,*) "I am all heater length"
        else 
            write(00,*) "Error: Please check circularGeo and heaterLength!"
            call MPI_finalize(ierr)
            stop
        endif
        write(00,*) 'Diameter=',real(2.0d0*radius0)
        if(dabs(lengthUnit-diameter0).GT.1e-6) then
            write(00,*) "WARNING: lengthUnit should be equal to diameter0 for circularGeo"
            write(00,*) "Error: Please check the geometry setting of the cell!"
            call MPI_finalize(ierr)
            stop
        endif
    elseif(flowGeometry.EQ.cornerLessGeo) then
        write(00,*) "I am cornerLessGeo"
        write(00,*) "gapX=",gapX,"; gapY=",gapY,"; thickness=",thickness
    elseif(flowGeometry.EQ.squareGeo) then
        write(00,*) "WARNING: flowGeometry should be irregularGeo; while its actual value is ", flowGeometry
        write(00,*) "Error: Please check the geometry setting of the cell!"
        call MPI_finalize(ierr)
        stop
    endif
#endif
#ifndef irregularGeo
    write(00,*) "I am regular geometry (Square Cavity)"            
    if(flowGeometry.NE.squareGeo) then
        write(00,*) "WARNING: flowGeometry should be", squareGeo, "; while its actual value is ", flowGeometry
        write(00,*) "Error: Please check the geometry setting of the cell!"
        call MPI_finalize(ierr)
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
        call MPI_finalize(ierr)
        stop
    endif

close(00)
endif
    
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
    
    obst = 0
    rho = 1.0d0
    
#ifdef irregularGeo
    if(myID.EQ.0) then
        if(flowGeometry.EQ.circularGeo) then
            obstGlobal = 1
            fluidNumMax = 0
            do j=1,ny
                do i=1,nx
                    if( ((dble(i)-xCenter)**2.0d0+(dble(j)-yCenter)**2.0d0).LE.radius0**2.0d0 ) then
                        obstGlobal(i,j) = 0
                        fluidNumMax = fluidNumMax+1
                    endif
                enddo
            enddo
        elseif(flowGeometry.EQ.cornerLessGeo) then
            obstGlobal = 0
            fluidNumMax = nx*ny
            do j=1+gapY,ny-gapY
                do i=1+gapX,nxFourth+gapX
                    if( ((i+j).GE.(nyFourth+gapX+gapY)).AND.((i+j).LE.(nyFourth+gapX+gapY+thickness)) ) then
                        obstGlobal(i,j) = 1 !~!~Bottom left corner
                        fluidNumMax = fluidNumMax-1
                    endif
                    if( ((-i+j).GE.(nyThreeFourths-gapX-gapY-thickness)).AND.((-i+j).LE.(nyThreeFourths-gapX-gapY)) ) then
                        obstGlobal(i,j) = 1 !~!~ Top left corner
                        fluidNumMax = fluidNumMax-1
                    endif
                enddo
            enddo
            do j=1+gapY,ny-gapY
                do i=nx-nxFourth+1-gapX,nx-gapX
                    if( ((-i+j).GE.(-nyThreeFourths+gapX+gapY)).AND.((-i+j).LE.(-nyThreeFourths+gapX+gapY+thickness)) ) then
                        obstGlobal(i,j) = 1 !~!~Bottom right corner
                        fluidNumMax = fluidNumMax-1
                    endif
                    if( ((i+j).GE.(ny+nyThreeFourths-gapX-gapY-thickness+1)).AND.((i+j).LE.(ny+nyThreeFourths-gapX-gapY+1)) ) then
                        obstGlobal(i,j) = 1 !~!~Top right corner
                        fluidNumMax = fluidNumMax-1
                    endif
                enddo
            enddo
        endif
    endif
    
    iRoot = 0
    call MPI_SCATTERV(obstGlobal, count2d, displ2d, MPI_integer4,  &
                                obstFit(1,1), count2d(myID), MPI_integer4, iRoot, MPI_COMM_WORLD, ierr)
    
    do j=0,nyLocal+1
        do i=1,nx
            obst(i,j) = obstFit(i,j)
        enddo
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !-------------------boundary data exchange
    call MPI_SENDRECV(obst(0,iEnd), (nx+2), MPI_integer4, rightNeighbor, 110, &
                                    obst(0,iStartMinus1), (nx+2), MPI_integer4, leftNeighbor, 110, &
                                    MPI_COMM_WORLD, iStatus, ierr)
    call MPI_SENDRECV(obst(0,iStart), (nx+2), MPI_integer4, leftNeighbor, 120, &
                                    obst(0,iEndPlus1), (nx+2), MPI_integer4, rightNeighbor, 120, &
                                    MPI_COMM_WORLD, iStatus, ierr) 
                                    
    deallocate (obstFit)
#endif

    if(loadInitField.EQ.0) then
        
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0
            
        if(myID.EQ.0) then
            open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
            write(00,*) "Initial field is set exactly"
            if(reloadDimensionlessTime.NE.0) then
                write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
                call MPI_finalize(ierr)
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
            if(flowGeometry.EQ.circularGeo) write(00,*) "Temperature B.C. for top/bottom curved walls are:===Hot/cold wall==="
            if(flowGeometry.EQ.cornerLessGeo) write(00,*)"Temperature B.C. for irregular curved walls are:===Adiabatic wall==="
#endif
        endif !~!~endif myID.EQ.0
        
#ifdef irregularGeo
        if(flowGeometry.EQ.circularGeo) then
            tempJ = 0
            do i=0,myID-1,1
                tempJ = tempJ+count1d(i)
            enddo
            do j=iStart,iEnd 
                do i=1,nx 
                    if(obst(i,j).EQ.0) then
                        T(i,j) = 0.0d0
                        do alpha=1,8
                            if(obst(i+ex(alpha),j+ey(alpha)).EQ.1) then
                            
                                if(heaterLength.EQ.oneFourthHeater) then
                                    if((j+tempJ).LE.int(yCenter-dsqrt(2.0d0)/2.0d0*radius0)) then
                                        T(i,j) = Thot
                                    endif
                                    if((j+tempJ).GE.int(yCenter+dsqrt(2.0d0)/2.0d0*radius0)) then
                                        T(i,j) = Tcold
                                    endif
                                elseif(heaterLength.EQ.oneThirdHeater) then
                                    if((j+tempJ).LE.int(yCenter-0.5d0*radius0)) then
                                        T(i,j) = Thot
                                    endif
                                    if((j+tempJ).GE.int(yCenter+0.5d0*radius0)) then
                                        T(i,j) = Tcold
                                    endif
                                elseif(heaterLength.EQ.allHeater) then
                                    T(i,j) = Thot-dble(j+tempJ)*dble(Thot-Tcold)/diameter0
                                endif
                                    
                            endif
                            
                        enddo
                    endif
                enddo
            enddo
        endif
#endif

#ifdef VerticalWallsConstT
        if(myID.EQ.0) then
            open(unit=00, file="SimulationSettings.txt", status='unknown', position='append')
            write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
            close(00)
        endif
        
        do j=iStart, iEnd
            !--Just walls--
            T(1,j) = Thot
            T(nx,j)= Tcold
            
            !~ !--Walls & bulk: linear profile
            !~ do i=1,nx
                !~ T(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
            !~ enddo
        enddo
#endif
    
#ifdef HorizontalWallsConstT
        if(myID.EQ.0) then
            open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
            write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
            close(00)
        endif
        
        do i=1,nx
            !~T(i,1) = Thot
            !~T(i,ny) = Tcold
            do j = iStart, iEnd
                T(i,j) = dble(j-1+start1d(myID)-1)/dble(ny-1)*(Tcold-Thot)+Thot
            enddo
        enddo
#endif

        if(myID.EQ.0) then
            open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
#ifdef VerticalWallsAdiabatic
            write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif
#ifdef HorizontalWallsAdiabatic
            write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif
            close(00)
        endif !~!~endif myID.EQ.0
    
        f = 0.0d0
        g = 0.0d0
        !$omp parallel do default(none) shared(f,g,u,v,rho,T,ex,ey,omega,omegaT,iStart,iEnd) private(i,j,alpha,us2,un)
        do j=iStart,iEnd
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
    
        if(myID.EQ.0) then
            open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
            if(reloadDimensionlessTime.EQ.0) then
                write(00,*) "Error: since loadInitField.EQ.1, reloadDimensionlessTime should not be 0"
                call MPI_finalize(ierr)
                stop
            endif
            write(00,*) "Load initial field from previous simulation: ../reloadFile/backupFile- >>>"
            write(reloadFileName, *) reloadbinFileNum
            reloadFileName = adjustl(reloadFileName)
            open(unit=01,file="../reloadFile/backupFile-"//trim(reloadFileName)//".bin",form="unformatted",access="sequential",status='old')
            read(01) (((fGlobal(alpha,i,j),alpha=0,8), i=1,nx), j=1,ny)
            read(01) (((gGlobal(alpha,i,j),alpha=0,4), i=1,nx), j=1,ny)
            read(01) ((uGlobal(i,j),i=1,nx), j=1,ny)
            read(01) ((vGlobal(i,j),i=1,nx), j=1,ny)
            read(01) ((tGlobal(i,j),i=1,nx), j=1,ny)
            close(01)
        
            dev = 0.0d0
            do j=1,ny
                do i=1,nx
                    rhoGlobal(i,j) = fGlobal(0,i,j)+fGlobal(1,i,j)+fGlobal(2,i,j)+fGlobal(3,i,j)+fGlobal(4,i,j)+fGlobal(5,i,j)+fGlobal(6,i,j)+fGlobal(7,i,j)+fGlobal(8,i,j)
                    dev = dev+gGlobal(0,i,j)+gGlobal(1,i,j)+gGlobal(2,i,j)+gGlobal(3,i,j)+gGlobal(4,i,j)-tGlobal(i,j)
                enddo
            enddo
            write(00,*)  "RELOAD: Deviation in temperature: ", real(dev)
            if(dev.GT.1.0d0) then
                write(00,*) "Error: too large Deviation when reload data!"
                call MPI_finalize(ierr)
                stop
            endif
            write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"
            close(00)
        endif !~!~end if myID.EQ.0

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
        iRoot = 0
        call MPI_SCATTERV(uGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION,  &
                                        u(1,1), count2d(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
                                        
        iRoot = 0
        call MPI_SCATTERV(vGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION,  &
                                        v(1,1), count2d(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
                                        
        iRoot = 0
        call MPI_SCATTERV(tGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION,  &
                                        T(1,1), count2d(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
                                        
        iRoot = 0
        call MPI_SCATTERV(rhoGlobal, count2d, displ2d, MPI_DOUBLE_PRECISION,  &
                                        rho(1,1), count2d(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
                                        
        iRoot = 0
        call MPI_SCATTERV(fGlobal, countD2Q9, displD2Q9, MPI_DOUBLE_PRECISION,  &
                                        f(0,1,1), countD2Q9(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
                                        
        iRoot = 0
        call MPI_SCATTERV(gGlobal, countD2Q5, displD2Q5, MPI_DOUBLE_PRECISION,  &
                                        g(0,1,1), countD2Q5(myID), MPI_DOUBLE_PRECISION, iRoot, MPI_COMM_WORLD, ierr)
    endif !~!~end if loadInitField 0 or 1
    
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

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,Fx,Fy,T,iStart,iEnd,obst,omega) private(i,j,alpha,s,m,m_post,meq,fSource) 
    do j=iStart,iEnd
        do i=1,nx
        
#ifdef irregularGeo
            if(obst(i,j).EQ.0) then
#endif

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

    f_post(0,i,j) = m_post(0)/9.0d0-m_post(1)/9.0d0+m_post(2)/9.0d0
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
                    
#ifdef irregularGeo
            else
                do alpha=0,8
                    f(alpha,i,j) = rho(i,j)*omega(alpha)
                    f_post(alpha,i,j) = rho(i,j)*omega(alpha)
                enddo
            endif
#endif

        enddo
    enddo
    !$omp end parallel do
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !-------------------boundary data exchange
call MPI_SENDRECV(f_post(0,-1,iEnd), 9*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 110, &
                                f_post(0,-1,iStartMinus1), 9*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 110, &
                                MPI_COMM_WORLD, iStatus, ierr)
call MPI_SENDRECV(f_post(0,-1,iStart), 9*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 120, &
                                f_post(0,-1,iEndPlus1), 9*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 120, &
                                MPI_COMM_WORLD, iStatus, ierr) 
                                
call MPI_SENDRECV(f_post(0,-1,iEnd-1), 9*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 130, &
                                f_post(0,-1,iStartMinus1-1), 9*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 130, &
                                MPI_COMM_WORLD, iStatus, ierr)
call MPI_SENDRECV(f_post(0,-1,iStart+1), 9*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 140, &
                                f_post(0,-1,iEndPlus1+1), 9*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 140, &
                                MPI_COMM_WORLD, iStatus, ierr)
    
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
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey,iStart,iEnd,obst) private(i,j,ip,jp,alpha)
    do j=iStart,iEnd
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
!*!*end subroutine streaming


!*!*begin subroutine calQ
    subroutine calQ(i,j,alpha,q)
    use commondata
    implicit none
    real(kind=8) :: i, j
    real(kind=8) :: q
    real(kind=8) :: x0, y0
    integer(kind=4) :: cNum
    integer(kind=4) :: alpha
    real(kind=8), parameter :: epsRadius=1e-9
    real(kind=8) :: qTemp

    q = 0.5d0
    qTemp = 0.5d0
    x0 = i+qTemp*dble(ex(alpha))
    y0 = j+qTemp*dble(ey(alpha))
    do while( dabs(dsqrt( (x0-dble(nxHalf))**2.0d0+(y0-dble(nyHalf))**2.0d0 )-radius0).GE.epsRadius) 
        if(dsqrt((x0-dble(nxHalf))**2.0d0+(y0-dble(nyHalf))**2.0d0).GT.radius0) then
            qTemp = qTemp/2.0d0
            x0 = x0-qTemp*dble(ex(alpha))
            y0 = y0-qTemp*dble(ey(alpha))
            q = q-qTemp
        elseif(dsqrt((x0-dble(nxHalf))**2.0d0+(y0-dble(nyHalf))**2.0d0).LT.radius0) then
            qTemp = qTemp/2.0d0
            x0 = x0+qTemp*dble(ex(alpha))
            y0 = y0+qTemp*dble(ey(alpha))
            q = q+qTemp
        else
            write(*,*) "error calQ!"
            call MPI_finalize(ierr)
            stop
        endif
    enddo

    return
    end subroutine calQ
!*!*begin subroutine calQ


!*!*Begin subroutine bounceback
    subroutine bounceback()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0, y0, q
    integer(kind=4) :: tempJ
    
#ifdef irregularGeo
    tempJ = 0
    do i=0,myID-1,1
        tempJ = tempJ+count1d(i)
    enddo
    
    !$omp parallel do default(none) shared(obst,f,f_post,ex,ey,r,istart,iEnd,tempJ) private(i,j,alpha,x0,y0,q)
    do j=iStart,iEnd
        do i=1,nx
            if(obst(i,j).EQ.0) then
            
                do alpha=1,8
                
                    if(obst(i+ex(alpha),j+ey(alpha)).EQ.1) then

                        if(flowGeometry.EQ.circularGeo) then
                            call calQ(dble(i),dble(j+tempJ),alpha,q)

                            if(q.LT.0.5d0) then
                                f(r(alpha),i,j) = q*(2.0d0*q+1.0d0)*f_post(alpha,i,j) &
                                                 +(1.0d0-4.0d0*q*q)*f_post(alpha,i-ex(alpha),j-ey(alpha)) &
                                                 -q*(1.0d0-2.0d0*q)*f_post(alpha,i-2*ex(alpha),j-2*ey(alpha))                               
                            elseif(q.GE.0.5d0) then
                                f(r(alpha),i,j) = f_post(alpha,i,j)/q/(2.0d0*q+1.0d0) &
                                                 +(2.0d0*q-1.0d0)/q*f_post(r(alpha),i,j) &
                                                 +(1.0d0-2.0d0*q)/(1.0d0+2.0d0*q) *f_post(r(alpha),i-ex(alpha),j-ey(alpha))
                            endif
                        
                        elseif(flowGeometry.EQ.cornerLessGeo) then
                            
                            f(r(alpha),i,j) = f_post(alpha,i,j)
                            
                        endif
                        
                    endif
    
                enddo
                
            endif
        enddo
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsPeriodicalU
    !$omp parallel do default(none) shared(f,f_post,iStart,iEnd) private(j)
    do j=iStart,iEnd
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
#endif

#ifdef VerticalWallsNoslip
    !$omp parallel do default(none) shared(f,f_post,iStart,iEnd) private(j)
    do j=iStart,iEnd
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
#endif

#ifdef HorizontalWallsNoslip
    if(myID.EQ.0) then
        !$omp parallel do default(none) shared(f,f_post,iStart) private(i,j)
        do i=1,nx
            !Bottom side (j=1)
            j = iStart
            f(2,i,j) = f_post(4,i,j)
            f(5,i,j) = f_post(7,i,j)
            f(6,i,j) = f_post(8,i,j)
        enddo
        !$omp end parallel do
    endif

    if(myID.EQ.nProc-1) then
        !$omp parallel do default(none) shared(f,f_post,iEnd) private(i,j)
        do i=1,nx
            !Top side (j=ny)
            j = iEnd
            f(4,i,j) = f_post(2,i,j)
            f(7,i,j) = f_post(5,i,j)
            f(8,i,j) = f_post(6,i,j)
        enddo
        !$omp end parallel do
    endif
#endif

    return
    end subroutine bounceback
!*!*End subroutine bounceback


!*!*Begin subroutine macro
    subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    
    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy,iStart,iEnd,obst) private(i,j)
    do j=iStart,iEnd
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

    !$omp parallel do default(none) shared(g,g_post,u,v,T,iStart,iEnd,obst,omegaT) private(i,j,alpha,n,neq,q,n_post) 
    do j=iStart,iEnd
        do i=1,nx
        
#ifdef irregularGeo
            if(obst(i,j).EQ.0) then
#endif
            
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
        
#ifdef irregularGeo
            else
                do alpha=0,4
                    g(alpha,i,j) = T(i,j)*omegaT(alpha)
                    g_post(alpha,i,j) = T(i,j)*omegaT(alpha)
                enddo
            endif
#endif

        enddo
    enddo
    !$omp end parallel do
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !-------------------boundary data exchange
call MPI_SENDRECV(g_post(0,-1,iEnd), 5*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 110, &
                                g_post(0,-1,iStartMinus1), 5*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 110, &
                                MPI_COMM_WORLD, iStatus, ierr)
call MPI_SENDRECV(g_post(0,-1,iStart), 5*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 120, &
                                g_post(0,-1,iEndPlus1), 5*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 120, &
                                MPI_COMM_WORLD, iStatus, ierr) 
                                
call MPI_SENDRECV(g_post(0,-1,iEnd-1), 5*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 130, &
                                g_post(0,-1,iStartMinus1-1), 5*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 130, &
                                MPI_COMM_WORLD, iStatus, ierr)
call MPI_SENDRECV(g_post(0,-1,iStart+1), 5*(nx+4), MPI_DOUBLE_PRECISION, leftNeighbor, 140, &
                                g_post(0,-1,iEndPlus1+1), 5*(nx+4), MPI_DOUBLE_PRECISION, rightNeighbor, 140, &
                                MPI_COMM_WORLD, iStatus, ierr) 
                                
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
    
    !$omp parallel do default(none) shared(g,g_post,ex,ey,iStart,iEnd,obst) private(i,j,ip,jp,alpha)
    do j=iStart,iEnd
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=0,4
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
    real(8) :: x0,y0,q,Cd1,Cd2,Cd3,Cd4
    integer(kind=4) :: tempJ
    
#ifdef irregularGeo
    tempJ = 0
    do i=0,myID-1,1
        tempJ = tempJ+count1d(i)
    enddo
    
    !$omp parallel do default(none) shared(obst,T,g,g_post,ex,ey,r,iStart,iEnd,tempJ) private(i,j,alpha,x0,y0,q,Cd1,Cd2,Cd3,Cd4)
    do j=iStart,iEnd
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=1,4

                    if(obst(i+ex(alpha),j+ey(alpha)).EQ.1) then

                        if(flowGeometry.EQ.circularGeo) then
                        
                            call calQ(dble(i),dble(j+tempJ),alpha,q)
                            
                            Cd1 = -1.0d0
                            Cd2 = (2.0d0*q-1.0d0)/(2.0d0*q+1.0d0)
                            Cd3 = (2.0d0*q-1.0d0)/(2.0d0*q+1.0d0)
                            Cd4 = 2.0d0/(2.0d0*q+1.0d0)
                            
                            if(heaterLength.EQ.oneFourthHeater) then
                
                                if((j+tempJ).LE.int(yCenter-dsqrt(2.0d0)/2.0d0*radius0)) then
                                    g(r(alpha),i,j) = Cd1*g_post(alpha,i,j)+Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) &
                                    +Cd4*(4.0d0+paraA)/10.0d0*Thot       
                                elseif((j+tempJ).GE.int(yCenter+dsqrt(2.0d0)/2.0d0*radius0)+1) then
                                    g(r(alpha),i,j) = Cd1*g_post(alpha,i,j)+Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) &
                                    +Cd4*(4.0d0+paraA)/10.0d0*Tcold       
                                else
                                    g(r(alpha),i,j) = -Cd1*g_post(alpha,i,j)-Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) 
                                endif
                               
                            elseif(heaterLength.EQ.oneThirdHeater) then

                                if((j+tempJ).LE.int(yCenter-0.5d0*radius0)) then
                                    g(r(alpha),i,j) = Cd1*g_post(alpha,i,j)+Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) &
                                    +Cd4*(4.0d0+paraA)/10.0d0*Thot       
                                elseif((j+tempJ).GE.int(yCenter+0.5d0*radius0)+1) then
                                    g(r(alpha),i,j) = Cd1*g_post(alpha,i,j)+Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) &
                                    +Cd4*(4.0d0+paraA)/10.0d0*Tcold       
                                else
                                    g(r(alpha),i,j) = -Cd1*g_post(alpha,i,j)-Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) 
                                endif

                            elseif(heaterLength.EQ.allHeater) then
                                    g(r(alpha),i,j) = Cd1*g_post(alpha,i,j)+Cd2*g_post(alpha,i-ex(alpha),j-ey(alpha))+Cd3*g_post(r(alpha),i,j) &
                                    +Cd4*(4.0d0+paraA)/10.0d0* &
                                    (Thot-dble(j+tempJ)*dble(Thot-Tcold)/diameter0)

                            endif
                        
                        elseif(flowGeometry.EQ.cornerLessGeo) then
                            
                                g(r(alpha),i,j) = g_post(alpha,i,j)
                        
                        endif
                   
                    endif
                   
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do
#endif

#ifdef HorizontalWallsAdiabatic
    if(myID.EQ.0) then
        !$omp parallel do default(none) shared(g,g_post,iStart) private(i,j)
        do i=1,nx
            !Bottom side (j=1)
            j = iStart
            g(2,i,j) = g_post(4,i,j)
        enddo
        !$omp end parallel do
    endif

    if(myID.EQ.nProc-1) then
        !$omp parallel do default(none) shared(g,g_post,iEnd) private(i,j)
        do i=1,nx
            !Top side (j=ny)
            j = iEnd
            g(4,i,j) = g_post(2,i,j)
        enddo
        !$omp end parallel do
    endif
#endif

#ifdef HorizontalWallsConstT
    if(myID.EQ.0) then
        !$omp parallel do default(none) shared(g,g_post,iStart) private(i,j)
        do i=1,nx
            !Bottom side (j=1)
            j = iStart
            g(2,i,j) = -g_post(4,i,j)+(4.0d0+paraA)/10.0d0*Thot
        enddo
        !$omp end parallel do
    endif

    if(myID.EQ.nProc-1) then
        !$omp parallel do default(none) shared(g,g_post,iEnd) private(i,j)
        do i=1,nx
            !Top side (j=ny)
            j = iEnd
            g(4,i,j) = -g_post(2,i,j)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
        !$omp end parallel do
    endif
#endif

#ifdef VerticalWallsConstT
    !$omp parallel do default(none) shared(g,g_post,iStart,iEnd) private(j) 
    do j=iStart,iEnd
        !Left side
        g(1,1,j) = -g_post(3,1,j)+(4.0d0+paraA)/10.0d0*Thot

        !Right side
        g(3,nx,j) = -g_post(1,nx,j)+(4.0d0+paraA)/10.0d0*Tcold
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsAdiabatic
    !$omp parallel do default(none) shared(g,g_post,iStart,iEnd) private(j) 
    do j=iStart,iEnd
        !Left side
        g(1,1,j) = g_post(3,1,j)

        !Right side
        g(3,nx,j) = g_post(1,nx,j)
    enddo
    !$omp end parallel do
#endif

#ifdef VerticalWallsPeriodicalT
    !$omp parallel do default(none) shared(g,g_post,iStart,iEnd) private(j) 
    do j=iStart,iEnd
        !Left side
        g(1,1,j) = g_post(1,nx,j)

        !Right side
        g(3,nx,j) = g_post(3,1,j)
    enddo
    !$omp end parallel do
#endif

    return
    end subroutine bouncebackT
!*!*end subroutine bouncebackT


!*!*Begin subroutine macroT
    subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    
    !$omp parallel do default(none) shared(g,T,iStart,iEnd,obst) private(i,j)
    do j=iStart,iEnd
        do i=1,nx
            if(obst(i,j).EQ.0) then
                T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            endif
        enddo
    enddo
    !$omp end parallel do

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
    real(kind=8) :: error1All, error2All, error5All, error6All
    
    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,T,Tp,iStart,iEnd,obst) private(i,j) reduction(+:error1,error2,error5,error6)
    do j=iStart,iEnd
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
    
    call MPI_ALLREDUCE(error1, error1All, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(error2, error2All, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(error5, error5All, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(error6, error6All, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    errorU = dsqrt(error1All)/dsqrt(error2All)
    errorT = error5All/error6All

    if(myID.EQ.0) then
        open(unit=01,file='convergence.log',status='unknown',position='append')
        write(01,*) itc,' ',errorU,' ',errorT
        close(01)
    endif

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
    write(03) ((uGlobal(i,j),i=1,nx),j=1,ny)
    write(03) ((vGlobal(i,j),i=1,nx),j=1,ny)
    write(03) ((tGlobal(i,j),i=1,nx),j=1,ny)
    close(03)

    return
    end subroutine output_binary
!*!*end subroutine output_binary

    
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
    write(05) (((fGlobal(alpha,i,j), alpha=0,8), i=1,nx), j=1,ny)
    write(05) (((gGlobal(alpha,i,j), alpha=0,4), i=1,nx), j=1,ny)
    write(05) ((uGlobal(i,j), i=1,nx), j=1,ny)
    write(05) ((vGlobal(i,j), i=1,nx), j=1,ny)
    write(05) ((tGlobal(i,j), i=1,nx), j=1,ny)
    close(05)
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
    write(00,*) "Backup f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
    close(00)

    return
    end subroutine backupData
!*!*end subroutine backupData
    
    
!*!*Begin subroutine output_Tecplot
    subroutine output_Tecplot()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    character(len=100) :: filename

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

#ifdef steadyFlow
    write(filename,'(i12.12)') itc
#endif
#ifdef unsteadyFlow
    pltFileNum = pltFileNum+1
    write(filename,'(i12.12)') pltFileNum
#endif
    filename = adjustl(filename)

    open(41,file='buoyancyCavity-'//trim(filename)//'.plt', access='stream', form='unformatted')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--integer(kind=4) value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file
#ifdef irregularGeo
    write(41) 6
#endif
#ifndef irregularGeo
    write(41) 5
#endif

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
#ifdef irregularGeo
    V6='obst'
    call dumpstring(V6)
#endif

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
#ifdef irregularGeo
    write(41) 1
#endif

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
                write(41) real(uGlobal(i,j))
                write(41) real(vGlobal(i,j))
                write(41) real(tGlobal(i,j))
#ifdef irregularGeo
                write(41) real(obstGlobal(i,j))
#endif
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
    !$omp parallel do default(none) shared(uGlobal,vGlobal,obstGlobal) private(i,j) reduction(+:angularMomentum)
    do j=1,ny 
        do i=1,nx 
            if(obstGlobal(i,j).EQ.0) then
                angularMomentum = angularMomentum+(i-nxHalf)*vGlobal(i,j)-(j-nyHalf)*uGlobal(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do
    angularMomentum = angularMomentum/dble(fluidNumMax)
    
    open(unit=01,file="angularMomentum.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), angularMomentum
    close(01)
    
    NuVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(vGlobal,tGlobal,obstGlobal) private(i,j) reduction(+:NuVolAvg_temp)
    do j=1,ny
        do i=1,nx
            if(obstGlobal(i,j).EQ.0) then
                NuVolAvg_temp = NuVolAvg_temp+vGlobal(i,j)*tGlobal(i,j)
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
    !$omp parallel do default(none) shared(uGlobal,vGlobal,obstGlobal) private(i,j) reduction(+:ReVolAvg_temp)
    do j=1,ny
        do i=1,nx       
            if(obstGlobal(i,j).EQ.0) then
                ReVolAvg_temp = ReVolAvg_temp+(uGlobal(i,j)*uGlobal(i,j)+vGlobal(i,j)*vGlobal(i,j))
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

