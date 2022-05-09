subroutine collision()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)


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
    
    return
end subroutine collision




subroutine streaming()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
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