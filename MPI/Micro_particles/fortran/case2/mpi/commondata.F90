module commondata
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)

    integer, parameter :: nx=101,ny=801
    integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

    real(8), parameter :: l0=4.0d0/1000.0d0 !! unit: cm
    real(8), parameter :: t0=8.0d0/100000.0d0 !! unit: s
    real(8), parameter :: tMax=2.5d0 !!unit: s

    integer :: itc
    integer, parameter :: itc_max=INT(tMax/t0)
    integer, parameter :: cNumMax=1
    
    real(8) :: xCenter(cNumMax), yCenter(cNumMax)
    real(8) :: xCenterOld(cNumMax), yCenterOld(cNumMax) 
    
    real(8) :: rationalOmega(cNumMax)
    real(8) :: rationalOmegaOld(cNumMax)
    
    real(8) :: Uc(cNumMax), Vc(cNumMax)
    real(8) :: UcOld(cNumMax), VcOld(cNumMax)
    
    real(8), parameter :: radius0=25.0d0/2.0d0
    real(8) :: radius(cNumMax)
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: rhoSolid=1.03d0
    real(8) :: rhoAvg

    real(8), parameter :: viscosity=0.05d0
    real(8), parameter :: tauf=3.0d0*viscosity+0.5d0
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)

    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    real(8) :: X(nx), Y(ny)
    real(8) :: u(nx,ny), v(nx,ny), rho(nx,ny)
    real(8) :: up(nx,ny), vp(nx,ny)
    real(8) :: f(0:8,-1:nx+2,-1:ny+2), f_post(0:8,-1:nx+2,-1:ny+2)
    real(8) :: omega(0:8)
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
    integer :: obst(0:nx+1,0:ny+1)
    integer :: obstNew(0:nx+1,0:ny+1)
    integer :: r(1:8)
    data r/3, 4, 1, 2, 7, 8, 5, 6/
    
    real(8) :: wallTotalForceX(cNumMax), wallTotalForceY(cNumMax)
    real(8) :: totalTorque(cNumMax)

    real(8), parameter :: gravity=980.0d0*t0**2.0d0/l0
    integer :: shiftTimes
end module commondata