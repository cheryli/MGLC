!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


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

end module commondata