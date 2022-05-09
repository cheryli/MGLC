!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


!~!~For steadyFlow, we define array up, vp, Tp to check  convergence;
!~!~for unsteadyFlow, we did not define array up, vp, Tp to save memory.
#define steadyFlow    
! #define unsteadyFlow

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
! #define RayleighBenardCell
! #define HorizontalWallsConstT
! #define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT
!!!!~~temperature B.C.~~

!!!~~temperature B.C. (for Side Heated Cell)~~
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!!!~~temperature B.C.~~


module ioFolder
    character(len=100) :: binFolderPrefix="./binFile/buoyancyCavity"
    character(len=100) :: pltFolderPrefix="./pltFile/buoyancyCavity"
    character(len=100) :: porousGeoFile="./case1/obstSquareCylinder-126.bin"
    character(len=100) :: particlePositionFolderPrefix="./particlePositionFile/particlePosition"
    character(len=100) :: particleReynoldsFolderPrefix="./particleReynoldsFile/particleReynolds"
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
    
    integer(kind=4), parameter :: total_nx=201, total_ny=201     !----Section 1----!
    integer :: nx, ny
    integer(kind=4), parameter :: flowGeometry=1    !----Section 1----!
    real(kind=8), parameter :: lengthUnit=dble(total_ny)     !----Section 1----!
    
    real(kind=8), parameter :: Rayleigh=1e7        !----Section 2----!
    real(kind=8), parameter :: Prandtl=0.71d0       !----Section 2----!
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

    real(kind=8), parameter :: shearReynolds = 0.0d0      !----Section 4----!
    
    integer(kind=4), parameter :: outputBinFile=0, outputPltFile=1                 !----Section *----!
            
    real(kind=8), parameter :: Pi=4.0d0*datan(1.0d0)
    real(kind=8), parameter :: tiltedAngle=dble(0.0d0/180.0d0*Pi)
    !-----------------------------------------------------------------------------------------------        
    !----Section 1----!
    integer(kind=4), parameter :: nxHalf=(total_nx-1)/2+1, nyHalf=(total_ny-1)/2+1
    integer(kind=4), parameter :: nxFourth=(total_nx-1)/4+1, nyFourth=(total_ny-1)/4+1
    integer(kind=4), parameter :: nxThreeFourths=3*(total_nx-1)/4+1, nyThreeFourths=3*(total_ny-1)/4+1
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
    real(kind=8), parameter :: U0=shearReynolds*viscosity/dble(total_ny)
    real(kind=8), parameter :: UwallTopLeft=U0, UwallTopRight=-U0, UwallBottomLeft=-U0, UwallBottomRight=U0
    real(kind=8), parameter :: UwallLeftTop=U0, UwallLeftBottom=U0, UwallRightTop=U0, UwallRightBottom=U0  !----Section 4----!!----Section 4----!
    !-----------------------------------------------------------------------------------------------
    !----Section 5----!
    real(kind=8) :: xp(0:total_nx+1), yp(0:total_ny+1)
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
    !-----------------------------------------------------------------------------------------------       
    
    ! mpi data
    integer :: rc, rank, num_process
    integer :: dims(0:1) = (0, 0), coords(0:1)
    logical :: periods(0:1)
    data periods/2*.false./
    integer :: comm2d, rank2d
    integer :: f_row_x, f_column_y, g_row_x, g_column_y
    integer :: nbr_left, nbr_right, nbr_top, nbr_bottom
    integer :: cnr_top_left, cnr_top_right, cnr_bottom_left, cnr_bottom_right
    integer :: i_start_global, j_start_global
    integer :: f_req(30), g_req(8)
end module commondata
        
    

program main
    use mpi
    use omp_lib    
    use commondata
    implicit none
    real(kind=8) :: timeStart, timeEnd
    real(kind=8) :: timeStart2, timeEnd2
    real(8) :: start_time, end_time
    integer(kind=4) :: myMaxThreads
    INTEGER(kind=4) :: time
    character(len=24) :: ctime, string
    integer :: name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
        
    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    call MPI_Dims_create(num_process, 2, dims, rc)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, rc)
    if(rank == 0) then
        write(*,*) "dimens is x*y = ", dims(0), "x", dims(1)
    endif

    ! get my new rank in decomposition
    call MPI_Comm_rank(comm2d, rank2d, rc)
    ! write(*,*) "process ", rank2d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    call MPI_Cart_get(comm2d, 2, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    ! write(*,*) "coords = ", coords(1), coords(2)
    ! write(*,*) "nx*ny = ", nx, ny

    ! get the neighbors
    call MPI_Cart_shift(comm2d, 0, 1, nbr_left, nbr_right, rc)
    call MPI_Cart_shift(comm2d, 1, 1, nbr_bottom, nbr_top, rc)
    ! write(*,*) "I'm process ", rank2d, "My neighbor lrbt is", nbr_left, nbr_right, nbr_bottom, nbr_top
    call MPI_Cart_find_corners()


    ! construct the datatype for the exchange 
    ! row_x exchange in y direction(non-contiguous) -- top and bottom
    call MPI_Type_vector(nx, 1, 9, MPI_REAL8, f_row_x, rc)
    call MPI_Type_commit(f_row_x, rc)
    ! column_y exchange in x direction(non-contiguous) -- left and right
    call MPI_Type_vector(ny, 1, 9 * (nx+2), MPI_REAL8, f_column_y, rc)
    call MPI_Type_commit(f_column_y, rc)
    ! row_x exchange in y direction(non-contiguous) -- top and bottom
    call MPI_Type_vector(nx, 1, 5, MPI_REAL8, g_row_x, rc)
    call MPI_Type_commit(g_row_x, rc)
    ! column_y exchange in x direction(non-contiguous) -- left and right
    call MPI_Type_vector(ny, 1, 5 * (nx+2), MPI_REAL8, g_column_y, rc)
    call MPI_Type_commit(g_column_y, rc)


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

    ! call output()
                
    call CPU_TIME(timeStart)
#ifdef _OPENMP
    timeStart2 = OMP_get_wtime()
#endif

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max).AND.(statisticallyConvergeState.EQ.0) )

        itc = itc+1
        
        call f_collision_with_message_exchange()

        call g_collision_with_message_exchange()

        call streaming()

        call bounceback()

        call macro()

        call streamingT()

        call bouncebackT()
        
        call macroT()

#ifdef steadyFlow
        if(MOD(itc,2000).EQ.0) call check()
#endif

        ! timer test
        if (mod(itc, 10000) == 0) then 
            exit
        endif

! #ifdef SideHeatedCell
        ! if( (MOD(itc,backupInterval*int(timeUnit)).EQ.0).AND.(outputBinFile.EQ.1) ) call backupData()
! #endif


! #ifdef RayleighBenardCell
!         if( MOD(itc, int(outputFrequency*timeUnit)).EQ.0 ) then
            
! #ifdef steadyFlow
!             ! if( (outputPltFile.EQ.1).AND.(MOD(itc, backupInterval*int(timeUnit)).EQ.0) ) call output_Tecplot()
! #endif
            
!             ! if(statisticallyStationaryState.EQ.1) then
!             !     if(outputBinFile.EQ.1) then
!             !         ! call output_binary()
!             !         ! if(MOD(itc, int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit)).EQ.0) call backupData()
!             !     endif   
!             !     ! if(outputPltFile.EQ.1) then
!             !     !     call output_Tecplot()
!             !     ! endif
!             ! endif
            
!             call calNuRe()
            
!         endif
! #endif
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()

    if(rank == 0) then
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
    endif
    
    ! if(outputBinFile.EQ.1) call backupData()
    
    ! call output()
        
    call CPU_TIME(timeEnd)
#ifdef _OPENMP
    timeEnd2 = OMP_get_wtime()
#endif
    
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
#ifdef _OPENMP
    write(00,*) "Time (OMP) = ", real(timeEnd2-timeStart2), "s"
#endif
    write(00,*) "Time (MPI) = ", real(end_time - start_time), "s"
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
    write(00,*) "    "

    write(00,*) "Successfully: DNS completed!"
    
    string = ctime( time() )
    write(00,*) 'End:   ', string
    close(00)

    call MPI_Finalize(rc)

    
    stop
end program main



subroutine decompose_1d(total_n, local_n, rank, num_process, i_start_global)
    implicit none
    integer, intent(in) :: total_n, rank, num_process
    integer, intent(out) :: local_n, i_start_global

    local_n = total_n / num_process

    if (rank < MOD(total_n, num_process)) then
        local_n = local_n + 1
    endif

    if (local_n > total_n / num_process) then ! --- 5 5 '5' 4 4 4
        i_start_global = local_n * rank
    else                    ! --- 5 5 5 4 '4' 4
        i_start_global = local_n * rank + mod(total_n, num_process)
    endif

end subroutine decompose_1d


subroutine MPI_Cart_find_corners()
    use mpi
    use commondata
    implicit none

    call MPI_Cart_shift_2d(1, 1, cnr_top_right)
    call MPI_Cart_shift_2d(1, -1, cnr_bottom_right)
    call MPI_Cart_shift_2d(-1, 1, cnr_top_left)
    call MPI_Cart_shift_2d(-1, -1, cnr_bottom_left)

end subroutine MPI_Cart_find_corners


subroutine MPI_Cart_shift_2d(idx0, idx1, corner_rank)
    use mpi
    use commondata
    implicit none
    integer, intent(in) :: idx0, idx1
    integer, intent(out) :: corner_rank 
    integer :: new_coords(0:1)

    new_coords(0) = coords(0) + idx0
    new_coords(1) = coords(1) + idx1
    if (new_coords(0) < 0 .OR. new_coords(0) > dims(0)-1) then  
        ! beyond the left/right boundary
        if (periods(0) == .false.) then
            corner_rank = MPI_PROC_NULL
            return
        else
            new_coords(0) = mod(new_coords(0) + dims(0), dims(0))
        endif
    else if (new_coords(1) < 0 .OR. new_coords(1) > dims(1)-1) then
        ! beyond the top/bottom boundary
        if (periods(1) == .false.) then
            corner_rank = MPI_PROC_NULL
            return
        else
            new_coords(1) = mod(new_coords(1) + dims(1), dims(1))
        endif
    endif

    call MPI_Cart_rank(comm2d, new_coords, corner_rank, rc)

end subroutine MPI_Cart_shift_2d


 
subroutine initial()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j, start, end
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2
    character(len=100) :: reloadFileName
    real(kind=8) :: dev
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    fluidNumMax = total_nx*total_ny
    
if (rank == 0) then
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
    write(00,*) 'Mesh:',total_nx, total_ny
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
    if(dabs(lengthUnit-dble(total_ny)).GT.1e-6) then
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
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif
#ifdef HorizontalWallsConstT
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
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

close(00)

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

        ! do i=1,nxHalf
        !     u(i,ny) = UwallTopLeft
        !     u(i,1)  = UwallBottomLeft
        ! enddo

        ! if (i_start_global + nx  <= nxHalf) then
        !     end = nx
        ! elseif (i_start_global >= nxHalf) then
        !     end = 0
        ! else
        !     end = nxHalf - i_start_global
        ! endif
        ! do i=1,end
        !     u(i,ny) = UwallTopLeft
        !     u(i,1)  = UwallBottomLeft
        ! enddo

        ! do i=nxHalf+1, nx
        !     u(i,ny) = UwallTopRight
        !     u(i,1) = UwallBottomRight
        ! enddo
        ! do j=2,nyHalf
        !     v(1,j) = UwallLeftBottom
        !     v(nx,j) = UwallRightBottom
        ! enddo
        ! do j=nyHalf+1, ny-1
        !     v(1,j) = UwallLeftTop
        !     v(nx,j) = UwallRightTop
        ! enddo


#ifdef VerticalWallsConstT
    do j=1,ny
        !~T(1,j) = Thot
        !~T(nx,j) = Tcold
        do i=1,nx
            ! T(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
            T(i, j) = dble(i_start_global + i - 1) / dble(total_nx - 1) * (Tcold - Thot) + Thot
        enddo
    enddo
#endif
#ifdef HorizontalWallsConstT
    do i=1,nx
        !~T(i,1) = Thot
        !~T(i,ny) = Tcold
        do j=1,ny
            ! T(i,j) = dble(j-1)/dble(ny-1)*(Tcold-Thot)+Thot
            T(i, j) = dble(j_start_global + j - 1) / dble(total_ny - 1) * (Tcold - Thot) + Thot
        enddo
    enddo
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

    ! elseif(loadInitField.EQ.1) then
    !     if(reloadDimensionlessTime.EQ.0) then
    !         write(00,*) "Error: since loadInitField.EQ.1, reloadDimensionlessTime should not be 0"
    !         stop
    !     endif
    !     write(00,*) "Load initial field from previous simulation >>>"
    !     write(reloadFileName, *) reloadbinFileNum
    !     reloadFileName = adjustl(reloadFileName)
    !     open(unit=01,file="../reloadFile/backupFile-"//trim(reloadFileName)//".bin",form="unformatted",access="sequential",status='old')
    !     read(01) (((f(alpha,i,j), alpha=0,8), i=1,nx), j=1,ny)
    !     read(01) (((g(alpha,i,j), alpha=0,4), i=1,nx), j=1,ny)
    !     read(01) ((u(i,j),i=1,nx), j=1,ny)
    !     read(01) ((v(i,j),i=1,nx), j=1,ny)
    !     reaD(01) ((T(i,j),i=1,nx), j=1,ny)
    !     close(01)
        
    !     dev = 0.0d0
    !     do j=1,ny
    !         do i=1,nx
    !             rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    !             dev = dev+g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    !         enddo
    !     enddo
    !     write(00,*) "RELOAD: Deviation in temperature: ", real(dev)
    !     if(dev.GT.1.0d0) then
    !         write(00,*) "Error: too large Deviation when reload data!"
    !         stop
    !     endif
    !     write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"
    else
        if (rank == 0) then
            write(00,*) "Error: initial field is not properly set"
        endif
        stop
    endif
    

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

    
subroutine f_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none

    ! calculate the boundary first
    call collision(1, nx, 1, 1)     ! bottom voundary
    call collision(1, nx, ny, ny)   ! top boundary
    call collision(1, 1, 1, ny)     ! left boundary
    call collision(nx, nx, 1, ny)   ! right boundary

    ! start meassage exchange
    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(f_post(2, 1, ny), 1, f_row_x, nbr_top, 2, comm2d, f_req(1), rc)
    call MPI_Irecv(f_post(2, 1, 0), 1, f_row_x, nbr_bottom, 2, comm2d, f_req(2), rc)

    call MPI_Isend(f_post(5, 1, ny), 1, f_row_x, nbr_top, 5, comm2d, f_req(3), rc)
    call MPI_Irecv(f_post(5, 1, 0), 1, f_row_x, nbr_bottom, 5, comm2d, f_req(4), rc)

    call MPI_Isend(f_post(6, 1, ny), 1, f_row_x, nbr_top, 6, comm2d, f_req(5), rc)
    call MPI_Irecv(f_post(6, 1, 0), 1, f_row_x, nbr_bottom, 6, comm2d, f_req(6), rc)

    ! message passing to bottom(j--)
    call MPI_Isend(f_post(4, 1, 1), 1, f_row_x, nbr_bottom, 4, comm2d, f_req(7), rc)
    call MPI_Irecv(f_post(4, 1, ny+1), 1, f_row_x, nbr_top, 4, comm2d, f_req(8), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, f_row_x, nbr_bottom, 7, comm2d, f_req(9), rc)
    call MPI_Irecv(f_post(7, 1, ny+1), 1, f_row_x, nbr_top, 7, comm2d, f_req(10), rc)

    call MPI_Isend(f_post(8, 1, 1), 1, f_row_x, nbr_bottom, 8, comm2d, f_req(11), rc)
    call MPI_Irecv(f_post(8, 1, ny+1), 1, f_row_x, nbr_top, 8, comm2d, f_req(12), rc)


    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(f_post(1, nx, 1), 1, f_column_y, nbr_right, 1, comm2d, f_req(13), rc)
    call MPI_Irecv(f_post(1, 0, 1), 1, f_column_y, nbr_left, 1, comm2d, f_req(14), rc)

    call MPI_Isend(f_post(5, nx, 1), 1, f_column_y, nbr_right, 5, comm2d, f_req(15), rc)
    call MPI_Irecv(f_post(5, 0, 1), 1, f_column_y, nbr_left, 5, comm2d, f_req(16), rc)

    call MPI_Isend(f_post(8, nx, 1), 1, f_column_y, nbr_right, 8, comm2d, f_req(17), rc)
    call MPI_Irecv(f_post(8, 0, 1), 1, f_column_y, nbr_left, 8, comm2d, f_req(18), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(f_post(3, 1, 1), 1, f_column_y, nbr_left, 3, comm2d, f_req(19), rc)
    call MPI_Irecv(f_post(3, nx+1, 1), 1, f_column_y, nbr_right, 3, comm2d, f_req(20), rc)

    call MPI_Isend(f_post(6, 1, 1), 1, f_column_y, nbr_left, 6, comm2d, f_req(21), rc)
    call MPI_Irecv(f_post(6, nx+1, 1), 1, f_column_y, nbr_right, 6, comm2d, f_req(22), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, f_column_y, nbr_left, 7, comm2d, f_req(23), rc)
    call MPI_Irecv(f_post(7, nx+1, 1), 1, f_column_y, nbr_right, 7, comm2d, f_req(24), rc)   


    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Isend(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, comm2d, f_req(25), rc)
    call MPI_Irecv(f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, comm2d, f_req(26), rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Isend(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, comm2d, f_req(27), rc)
    call MPI_Irecv(f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, comm2d, f_req(28), rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Isend(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, comm2d, f_req(29), rc)
    call MPI_Irecv(f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, comm2d, f_req(30), rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Isend(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, comm2d, f_req(31), rc)
    call MPI_Irecv(f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, comm2d, f_req(32), rc)


    ! then calculate the inner points
    call collision(2, nx-1, 2, ny-1)

    ! ! wait untill collision_f finished
    ! call MPI_Waitall(32, f_req, MPI_STATUSES_IGNORE, rc)


end subroutine f_collision_with_message_exchange


subroutine g_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none

    ! calculate the boundary first
    call collisionT(1, nx, 1, 1)     ! bottom voundary
    call collisionT(1, nx, ny, ny)   ! top boundary
    call collisionT(1, 1, 1, ny)     ! left boundary
    call collisionT(nx, nx, 1, ny)   ! right boundary

    ! message tag:  0, 1, 2, discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(g_post(2, 1, ny), 1, g_row_x, nbr_top, 2, comm2d, g_req(1), rc)
    call MPI_Irecv(g_post(2, 1, 0), 1, g_row_x, nbr_bottom, 2, comm2d, g_req(2), rc)


   ! message passing to bottom(j--)
    call MPI_Isend(g_post(4, 1, 1), 1, g_row_x, nbr_bottom, 4, comm2d, g_req(3), rc)
    call MPI_Irecv(g_post(4, 1, ny+1), 1, g_row_x, nbr_top, 4, comm2d, g_req(4), rc)

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(g_post(1, nx, 1), 1, g_column_y, nbr_right, 1, comm2d, g_req(5), rc)
    call MPI_Irecv(g_post(1, 0, 1), 1, g_column_y, nbr_left, 1, comm2d, g_req(6), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(g_post(3, 1, 1), 1, g_column_y, nbr_left, 3, comm2d, g_req(7), rc)
    call MPI_Irecv(g_post(3, nx+1, 1), 1, g_column_y, nbr_right, 3, comm2d, g_req(8), rc)
                
    ! then calculate the inner points
    call collisionT(2, nx-1, 2, ny-1)

    ! ! wait untill collision_g finished
    ! call MPI_Waitall(8, g_req, MPI_STATUSES_IGNORE, rc)

end subroutine g_collision_with_message_exchange


subroutine collision(i_start, i_end, j_start, j_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)


    do j = j_start, j_end
        do i = i_start, i_end


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
    
    return
end subroutine collision


subroutine collisionT(i_start, i_end, j_start, j_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    do j = j_start, j_end
        do i = i_start, i_end
            
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




subroutine streaming()
    use mpi
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
    ! wait untill collision_f finished
    call MPI_Waitall(32, f_req, MPI_STATUSES_IGNORE, rc)

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
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0, y0, q
    

! #ifdef VerticalWallsPeriodicalU
!     !$omp parallel do default(none) shared(f, f_post) private(j)
!     do j=1,ny
!         !Left side (i=1)
!         f(1,1,j) = f_post(1,nx,j)
!         f(5,1,j) = f_post(5,nx,j)
!         f(8,1,j) = f_post(8,nx,j)
        
!         !Right side (i=nx)
!         f(3,nx,j) = f_post(3,1,j)
!         f(6,nx,j) = f_post(6,1,j)
!         f(7,nx,j) = f_post(7,1,j)
!     enddo
!     !$omp end parallel do
!     i = 1
!     j = 1
!     f(2, i, j) = f_post(2, nx, j)
!     f(6, i, j) = f_post(6, nx, j)
    
!     i = nx
!     j = 1
!     f(2, i, j) = f_post(2, 1, j)
!     f(5, i, j) = f_post(5, 1, j)
    
!     i = 1
!     j = ny
!     f(4, i, j) = f_post(4, nx, j)
!     f(7, i, j) = f_post(7, nx, j)
    
!     i = nx
!     j = ny
!     f(4, i, j) = f_post(4, 1, j)
!     f(8, i, j) = f_post(8, 1, j)
! #endif

! #ifdef VerticalWallsNoslip
!     !$omp parallel do default(none) shared(f, f_post, rho) private(j)
!     do j = 2, nyHalf
!         !Left side (i=1)
!         f(1,1,j) = f_post(3,1,j)
!         f(5,1,j) = f_post(7,1,j)-rho(1,j)*(-UwallLeftBottom)/6.0d0
!         f(8,1,j) = f_post(6,1,j)-rho(1,j)*UwallLeftBottom/6.0d0
        
!         !Right side (i=nx)
!         f(3,nx,j) = f_post(1,nx,j)
!         f(6,nx,j) = f_post(8,nx,j)-rho(nx,j)*(-UwallRightBottom)/6.0d0
!         f(7,nx,j) = f_post(5,nx,j)-rho(nx,j)*UwallRightBottom/6.0d0
!     enddo
!     !$omp end parallel do
!     !$omp parallel do default(none) shared(f, f_post, rho) private(j)
!     do j = nyHalf+1, ny-1
!         !Left side (i=1)
!         f(1,1,j) = f_post(3,1,j)
!         f(5,1,j) = f_post(7,1,j)-rho(1,j)*(-UwallLeftTop)/6.0d0
!         f(8,1,j) = f_post(6,1,j)-rho(1,j)*UwallLeftTop/6.0d0
        
!         !Right side (i=nx)
!         f(3,nx,j) = f_post(1,nx,j)
!         f(6,nx,j) = f_post(8,nx,j)-rho(nx,j)*(-UwallRightTop)/6.0d0
!         f(7,nx,j) = f_post(5,nx,j)-rho(nx,j)*UwallRightTop/6.0d0
!     enddo
!     !$omp end parallel do
! #endif

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

! #ifdef HorizontalWallsNoslip
!     !$omp parallel do default(none) shared(f, f_post, rho) private(i)
!     do i = 2, nxHalf
!         !Bottom side (j=1)
!         f(2,i,1) = f_post(4,i,1)
!         f(5,i,1) = f_post(7,i,1)-rho(i,1)*(-UwallBottomLeft)/6.0d0
!         f(6,i,1) = f_post(8,i,1)-rho(i,1)*UwallBottomLeft/6.0d0

!         !Top side (j=ny)
!         f(4,i,ny) = f_post(2,i,ny)
!         f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*UwallTopLeft/6.0d0
!         f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-UwallTopLeft)/6.0d0
!     enddo
!     !$omp end parallel do
!     !$omp parallel do default(none) shared(f, f_post, rho) private(i)
!     do i = nxHalf+1, nx-1
!         !Bottom side (j=1)
!         f(2,i,1) = f_post(4,i,1)
!         f(5,i,1) = f_post(7,i,1)-rho(i,1)*(-UwallBottomRight)/6.0d0
!         f(6,i,1) = f_post(8,i,1)-rho(i,1)*UwallBottomRight/6.0d0

!         !Top side (j=ny)
!         f(4,i,ny) = f_post(2,i,ny)
!         f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*UwallTopRight/6.0d0
!         f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-UwallTopRight)/6.0d0
!     enddo
!     !$omp end parallel do
! #endif

! #ifndef VerticalWallsPeriodicalU
!     !!~corners
!     i = 1
!     j = 1
!     f(1,i,j) = f_post(3,i,j)
!     f(2,i,j) = f_post(4,i,j)
!     f(6,i,j) = f_post(8,i,j)-rho(i,j)*UwallBottomLeft/6.0d0
!     f(8,i,j) = f_post(6,i,j)-rho(i,j)*UwallLeftBottom/6.0d0
!     f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallBottomLeft-UwallLeftBottom)/6.0d0
    
!     i = nx 
!     j = 1
!     f(2,i,j) = f_post(4,i,j)
!     f(3,i,j) = f_post(1,i,j)
!     f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallBottomRight)/6.0d0
!     f(7,i,j) = f_post(5,i,j)-rho(i,j)*UwallRightBottom/6.0d0
!     f(6,i,j) = f_post(8,i,j)-rho(i,j)*(UwallBottomRight-UwallRightBottom)/6.0d0

!     i = 1
!     j = ny 
!     f(1,i,j) = f_post(3,i,j)
!     f(4,i,j) = f_post(2,i,j)
!     f(5,i,j) = f_post(7,i,j)-rho(i,j)*(-UwallLeftTop)/6.0d0
!     f(7,i,j) = f_post(5,i,j)-rho(i,j)*UwallTopLeft/6.0d0
!     f(8,i,j) = f_post(6,i,j)-rho(i,j)*(-UwallTopLeft+UwallLeftTop)/6.0d0
    
!     i = nx 
!     j = ny 
!     f(3,i,j) = f_post(1,i,j)
!     f(4,i,j) = f_post(2,i,j)
!     f(6,i,j) = f_post(8,i,j)-rho(i,j)*(-UwallRightTop)/6.0d0
!     f(8,i,j) = f_post(6,i,j)-rho(i,j)*(-UwallTopRight)/6.0d0
!     f(7,i,j) = f_post(5,i,j)-rho(i,j)*(UwallTopRight+UwallRightTop)/6.0d0
!     !!~corners
! #endif

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


#ifdef VerticalWallsNoslip
    !Left side (i=1)
    if (coords(0) == 0) then
        do j=1,ny 
            f(1,1,j) = f_post(3,1,j)
            f(5,1,j) = f_post(7,1,j)
            f(8,1,j) = f_post(6,1,j)
        enddo
    endif

    !Right side (i=nx)
    if (coords(0) == dims(0)-1) then    
        do j=1,ny
            f(3,nx,j) = f_post(1,nx,j)
            f(6,nx,j) = f_post(8,nx,j)
            f(7,nx,j) = f_post(5,nx,j)
        enddo
    endif
#endif

#ifdef HorizontalWallsNoslip
    !Bottom side (j=1)
    if (coords(1) == 0) then
        do i=1,nx 
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)
            f(6,i,1) = f_post(8,i,1)
        enddo
    endif

    !Top side (j=ny)
    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            f(4,i,ny) = f_post(2,i,ny)
            f(7,i,ny) = f_post(5,i,ny)
            f(8,i,ny) = f_post(6,i,ny)
        enddo
    endif
#endif

    return
end subroutine bounceback

    

subroutine macro()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    do j=1,ny
        do i=1,nx
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*Fx(i,j) )/rho(i,j)
                v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*Fy(i,j) )/rho(i,j)
        enddo
    enddo

    return
end subroutine macro
    



subroutine streamingT()
    use mpi
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha


    ! wait untill collision_g finished
    call MPI_Waitall(8, g_req, MPI_STATUSES_IGNORE, rc)
    
    do j = 1, ny
        do i = 1, nx
                do alpha = 0, 4
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                
                    g(alpha,i,j) = g_post(alpha,ip,jp)
                enddo
        enddo
    enddo
    
    return
end subroutine streamingT


subroutine bouncebackT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: x0,y0,q,Cd1,Cd2,Cd3,Cd4


#ifdef HorizontalWallsAdiabatic
    !Bottom side
    if (coords(1) == 0) then
        do i = 1, nx 
            g(2, i, 1) = g_post(4, i, 1)
        enddo
    endif

    ! Top side
    if (coords(1) == dims(1) - 1) then
        do i = 1, nx 
            g(4, i, ny) = g_post(2, i, ny)
        enddo
    endif
#endif


#ifdef HorizontalWallsConstT
   !Bottom side
    if (coords(1) == 0) then
        do i = 1, nx 
            g(2, i, 1) = -g_post(4, i, 1)+(4.0d0+paraA)/10.0d0*Thot
        enddo
    endif

    ! Top side
    if (coords(1) == dims(1) - 1) then
        do i = 1, nx 
            g(4, i, ny) = -g_post(2, i, ny)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
    endif
#endif

#ifdef VerticalWallsAdiabatic
    !Left side
    if (coords(0) == 0) then
        do j = 1, ny 
            g(1,1,j) = g_post(3,1,j)
        enddo
    endif

    !Right side
    if (coords(0) == dims(0) - 1) then
        do j = 1, ny 
            g(3,nx,j) = g_post(1,nx,j)
        enddo
    endif
#endif

#ifdef VerticalWallsConstT
    !Left side
    if (coords(0) == 0) then
        do j = 1, ny 
            g(1,1,j) = -g_post(3,1,j)+(4.0d0+paraA)/10.0d0*Thot
        enddo
    endif

    !Right side
    if (coords(0) == dims(0) - 1) then
        do j = 1, ny 
            g(3,nx,j) = -g_post(1,nx,j)+(4.0d0+paraA)/10.0d0*Tcold
        enddo
    endif
#endif


    
! #ifdef VerticalWallsPeriodicalT
!     !$omp parallel do default(none) shared(g,g_post) private(j) 
!     do j=1,ny
!         !Left side
!         g(1,1,j) = g_post(1,nx,j)

!         !Right side
!         g(3,nx,j) = g_post(3,1,j)
!     enddo
!     !$omp end parallel do
! #endif

    return
end subroutine bouncebackT



subroutine macroT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: a
    
    do j = 1, ny
        do i = 1, nx
                T(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        enddo
    enddo

    return
end subroutine macroT
    

#ifdef steadyFlow
    subroutine check()
    use mpi
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    

    do j=1,ny
        do i=1,nx
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                error5 = error5+dABS( T(i,j)-Tp(i,j) )
                error6 = error6+dABS( T(i,j) )
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j)
                Tp(i,j) = T(i,j)
        enddo
    enddo

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        open(unit=01,file='convergence.log',status='unknown',position='append')
            write(01,*) itc,' ',errorU,' ',errorT
        close(01)
        
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    
    return
    end subroutine check
#endif


subroutine output()
    use mpi
    use commondata
    integer :: i, j
    integer :: p_rank, num(0:3) ,dx = 0, dy = 0, new_coords(0:1)
    real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :), total_T(:, :), stream(:, :), vorticity(:, :)
    real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :), tmp_T(:, :)

    if (rank2d > 0) then
        ! rank != 0 send data
        num(0) = nx
        num(1) = ny
        num(2) = i_start_global
        num(3) = j_start_global
        ! send to rank 0
        call MPI_Send(num, 4, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size and origion
        call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
        call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
        call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
        call MPI_Send(T, nx*ny, MPI_REAL8, 0, 4, comm2d, rc)
    else
        ! rank 0 collect data
        allocate(total_u(total_nx, total_ny))
        allocate(total_v(total_nx, total_ny))
        allocate(total_rho(total_nx, total_ny))
        allocate(total_T(total_nx, total_ny))

        dx = i_start_global
        dy = j_start_global

        ! collect data from rank 0
        do j = 1, ny
            do i = 1, nx
                total_u(dx + i, dy + j) = u(i, j)
                total_v(dx + i, dy + j) = v(i, j)
                total_rho(dx + i, dy + j) = rho(i, j)
                total_T(dx + i, dy + j) = T(i, j)
            enddo
        enddo

         ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) - 1
            ! receive the block size and origion
            call MPI_Recv(num, 4, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)
            ! creat buffer
            allocate(tmp_u(num(0), num(1)))
            allocate(tmp_v(num(0), num(1)))
            allocate(tmp_rho(num(0), num(1)))
            allocate(tmp_T(num(0), num(1)))
            
            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_T, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 4, comm2d, MPI_STATUS_IGNORE, rc)

            ! determine the origin
            dx = num(2)
            dy = num(3)

            ! assign data
            do j = 1, num(1)
                do i = 1, num(0)
                    total_u(dx + i, dy + j) = tmp_u(i, j)
                    total_v(dx + i, dy + j) = tmp_v(i, j)
                    total_rho(dx + i, dy + j) = tmp_rho(i, j)
                    total_T(dx + i, dy + j) = tmp_T(i, j)
                enddo
            enddo

            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_rho)
            deallocate(tmp_T)
        enddo


        allocate(stream(total_nx, total_ny))
        allocate(vorticity(total_nx, total_ny))

        call compute_stream_vorticity(stream, vorticity, total_u, total_v, total_nx, total_ny)
        call output_Tecplot(xp, yp, total_u, total_v, total_rho, total_T, stream, vorticity, total_nx, total_ny, itc)
        call output_binary(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        call out_Velocity_Nu(total_u, total_v, total_T, total_nx, total_ny, diffusivity, lengthUnit)

        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_rho)
        deallocate(total_T)
        deallocate(stream)
        deallocate(vorticity)

    endif


end subroutine output

subroutine out_Velocity_Nu(u, v, T, nx, ny, diffusivity, lengthUnit)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: diffusivity, lengthUnit
    real(8), intent(in) :: u(nx, ny), v(nx, ny), T(nx, ny)
    integer :: i, j, nx_half, ny_half
    real(8) :: max_u, max_v, u_i, v_i
    real(8) :: Nu, Nu_vol, Nu_0, Nu_half, Nu_max, Nu_min, Nu_max_id, Nu_min_id

    nx_half = (1 + nx) / 2
    ny_half = (1 + ny) / 2

    ! find maximun u along vertical mid-plane, and the maximun v along horizontal mid-plane
    
    max_v = -1e30
    do i = 1, nx
        if (v(i, ny_half) > max_v) then
            max_v = v(i, ny_half)
            v_i = i
        endif
    enddo

    max_u = -1e30
    do j = 1, ny
        if (u(nx_half, j) > max_u) then
            max_u = u(nx_half, j)
            u_i = j
        endif
    enddo

    max_u = max_u / (diffusivity / lengthUnit)
    max_v = max_v / (diffusivity / lengthUnit)
    u_i = u_i / ny
    v_i = v_i / nx

    open(unit=77,file="result_parameters.txt",status='unknown',position="append")
        write(77,300) max_u, u_i, max_v, v_i
        write(*,300) max_u, u_i, max_v, v_i
300     format('max_u = ', f13.7," at y = ", f13.7, 'max_v = ', f13.7, " at x = ", f13.7)
    close(77)

    

    ! calculate Nu 
    Nu_vol = 0.0d0
    Nu_0 = 0.0d0
    Nu_half = 0.0d0
    Nu_max = -1e30
    Nu_min = 1e30
    do j = 1, ny
        do i = 1, nx
            Nu = v(i, j) * T(i, j)
            Nu_vol = Nu_vol + Nu
            if (i == 1) then     !! hot wall
                Nu_0 = Nu_0 + Nu
                if (Nu > Nu_max) then
                    Nu_max = Nu
                    Nu_max_id = j
                endif
                if (Nu < Nu_min) then
                    Nu_min = Nu
                    Nu_min_id = j
                endif
            endif

            if (i == nx_half) then  !! vertical mid-plane
                Nu_half = Nu_half + Nu
            endif
        enddo
    enddo

    Nu_vol = Nu_vol / (nx * ny) * lengthUnit / diffusivity + 1.0
    Nu_0 = Nu_0 / ny * lengthUnit / diffusivity + 1.0
    Nu_half = Nu_half / ny * lengthUnit / diffusivity + 1.0
    Nu_max = Nu_max * lengthUnit / diffusivity + 1.0
    Nu_min = Nu_min * lengthUnit / diffusivity + 1.0
    Nu_max_id = Nu_max_id / ny
    Nu_min_id = Nu_min_id / ny


    open(unit=77,file="result_parameters.txt",status='unknown',position="append")
        write(77,301) Nu_vol, Nu_0, Nu_half
        write(77,302) Nu_max, Nu_max_id, Nu_min
        write(*,301) Nu_vol, Nu_0, Nu_half
        write(*,302) Nu_max, Nu_max_id, Nu_min, Nu_min_id
301     format('Nu_vol = ', f13.7, 'Nu_0 = ', f13.7, 'Nu_half = ', f13.7)
302     format('Nu_max = ', f13.7," at y = ", f13.7, 'Nu_min = ', f13.7, " at y = ", f13.7)
    close(77)


end subroutine out_Velocity_Nu


subroutine output_binary(u, v, rho, T, nx, ny, itc)
    use ioFolder
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny)
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



subroutine output_Tecplot(xp, yp, u, v, rho, T, stream, vorticity, nx, ny, itc)
    use ioFolder
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny), T(nx, ny), stream(nx, ny), vorticity(nx, ny)
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7,V8
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

    write(41) 8

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
    V7='stream'
    call dumpstring(V7)
    V8='vorticity'
    call dumpstring(V8)


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
                write(41) real(stream(i,j))
                write(41) real(vorticity(i,j))
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



! subroutine backupData()
!     use commondata
!     implicit none
!     integer(kind=4) :: i, j, alpha
!     character(len=100) :: filename

! #ifdef steadyFlow
!     write(filename,*) itc
! #endif
! #ifdef unsteadyFlow
!     if(loadInitField.EQ.0) write(filename,*) binFileNum
!     if(loadInitField.EQ.1) write(filename,*) binFileNum+reloadbinFileNum
! #endif
!     filename = adjustl(filename)

!     open(unit=05,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
!     write(05) (((f(alpha,i,j),alpha=0,8), i=1,nx), j=1,ny)
!     write(05) (((g(alpha,i,j),alpha=0,4), i=1,nx), j=1,ny)
!     write(05) ((u(i,j), i=1,nx), j=1,ny)
!     write(05) ((v(i,j), i=1,nx), j=1,ny)
!     write(05) ((T(i,j), i=1,nx), j=1,ny)
!     close(05)
    
!     open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')
!     write(00,*) "Backup  f, g, u, v, T to the file: backupFile-", trim(filename),".bin"
!     close(00)

!     return
! end subroutine backupData


subroutine compute_stream_vorticity(stream, vorticity, u, v, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: u(nx, ny), v(nx, ny)
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
        vorticity(i, ny) = -2.0d0 * u(i, ny-1)
    enddo

    vorticity(1, ny) = 0
    vorticity(nx, ny) = 0

    vorticity(1, 1) = 0
    vorticity(nx, 1) = 0

end subroutine compute_stream_vorticity