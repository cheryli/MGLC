module commondata
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)

    integer :: nx, ny
    integer, parameter :: total_nx=201, total_ny=801

    real(8), parameter :: l0=1.0d0/100.0d0 !! unit: cm
    real(8), parameter :: t0=5.0d0/10000.0d0 !! unit: s
    real(8), parameter :: m0=1.0d0/1000000.0d0 !! unit: g
    real(8), parameter :: tMax=7.0d0 !! unit: s

    integer :: itc
    integer, parameter :: itc_max=INT(tMax/t0)
    integer, parameter :: cNumMax=64
    
    real(8) :: xCenter(cNumMax), yCenter(cNumMax)
    real(8) :: xCenterOld(cNumMax), yCenterOld(cNumMax) 
    
    real(8) :: rationalOmega(cNumMax)
    real(8) :: rationalOmegaOld(cNumMax)
    
    real(8) :: Uc(cNumMax), Vc(cNumMax)
    real(8) :: UcOld(cNumMax), VcOld(cNumMax)
    
    real(8), parameter :: radius0=20.0d0/2.0d0
    real(8) :: radius(cNumMax)
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: rhoSolid=1.01d0
    real(8) :: rhoAvg

    real(8), parameter :: viscosity=0.05d0
    real(8), parameter :: tauf=3.0d0*viscosity+0.5d0
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)

    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    real(8), allocatable :: X(:), Y(:)
    real(8), allocatable :: u(:, :), v(:, :), rho(:, :)
    real(8), allocatable :: up(:, :), vp(:, :)
    real(8), allocatable :: f(:, :, :), f_post(:, :, :)
    integer, allocatable :: obst(:, :)
    integer, allocatable :: obstNew(:, :)

    real(8) :: omega(0:8)
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
    integer :: r(1:8)
    data r/3, 4, 1, 2, 7, 8, 5, 6/
    
    real(8) :: wallTotalForceX(cNumMax), wallTotalForceY(cNumMax)
    real(8) :: totalTorque(cNumMax)
    real(8), parameter :: gravity=980.0d0*t0**2.0d0/l0

    real(8), parameter :: thresholdWall=6.0d0
    real(8), parameter :: stiffWall=0.02d0
    
    real(8), parameter :: thresholdParticle=6.0d0
    real(8), parameter :: stiffParticle=0.08d0


    ! real(8), parameter :: thresholdWall=4.0d0
    ! real(8), parameter :: stiffWall=0.10d0
    
    ! real(8), parameter :: thresholdParticle=4.0d0
    ! real(8), parameter :: stiffParticle=0.10d0

    ! mpi data
    integer :: rc, rank, num_process
    integer :: dims(0:1) = (0, 0), coords(0:1)
    logical :: periods(0:1)
    data periods/.false., .false./
    integer :: comm2d, rank2d
    integer :: fp_row_x, fp_column_y, type_f
    integer :: type_f_x, type_fi_x, type_f_y, type_fi_y, type_fp_y, type_fpi_y
    integer :: type_square3_f, type_square2_f, type_square2_fp, type_square2_fpi
    integer :: nbr_left, nbr_right, nbr_top, nbr_bottom
    integer :: cnr_top_left, cnr_top_right, cnr_bottom_left, cnr_bottom_right
    integer :: i_start_global, j_start_global

    ! particle mask
    integer :: local_mask(cNumMax)

end module commondata