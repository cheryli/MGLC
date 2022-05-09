!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


module commondata
    implicit none
    
    integer, parameter :: loadInitField=0
    integer, parameter :: flowReversal=0  !! value 0: do NOT simulate reversal; value 1: do simulate reversal

    integer, parameter :: total_nx = 51, total_ny = total_nx, total_nz = total_nx
    integer :: nx, ny, nz
    integer, parameter :: nxHalf = (total_nx-1)/2+1, nyHalf = (total_ny-1)/2+1, nzHalf = (total_nz-1)/2+1
    
    real(kind=8), parameter :: Rayleigh=1e6
    real(kind=8), parameter :: Prandtl=0.71d0
    real(kind=8), parameter :: Mach=0.1d0

    real(kind=8), parameter :: tauf=0.5d0+Mach*dble(total_nz)*DSQRT(3.0d0*Prandtl/Rayleigh)
    real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
    real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
    real(kind=8), parameter :: diffusivity=viscosity/Prandtl
    
    real(kind=8), parameter :: Ekman=0.001d0
    real(kind=8), parameter :: omegaRatating=viscosity/2.0d0/Ekman/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0
    real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/dble(total_nz)
    real(kind=8), parameter :: gBeta=gBeta1/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: timeUnit=dsqrt(dble(total_nz)/gBeta)  !!dble(ny*ny)/diffusivity
    integer, parameter :: dimensionlessTimeMax=1000
    integer, parameter :: flowReversalTime=20000
    integer :: itc
    integer, parameter :: itc_max=INT(dimensionlessTimeMax*timeUnit)
    
    real(kind=8), parameter :: epsU=1e-8
    real(kind=8), parameter :: epsT=1e-8
    real(kind=8) :: errorU, errorT
    
    real(kind=8) :: xp(0:total_nx+1), yp(0:total_ny+1), zp(0:total_nz+1)
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
    

    ! mpi data
    integer :: rc, rank, num_process 
    integer :: dims(0:2) = (/ 0, 0, 0 /), coords(0:2)
    logical :: periods(0:2) = (/ .false., .false., .false. /)
    integer :: comm3d, rank3d
    integer :: nbr_surface(1:6)
    integer :: nbr_line(7:18) 
    integer :: f_surface_x, f_surface_y, f_surface_z  ! surface data perpendicular to axis x, y, z
    integer :: f_line_x, f_line_y, f_line_z       ! line data parallel to axis x, y, z
    integer :: g_surface_x, g_surface_y, g_surface_z  
    integer :: g_line_x, g_line_y, g_line_z
    integer :: i_start_global, j_start_global, k_start_global
 

end module commondata