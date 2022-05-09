subroutine collisionT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

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


subroutine streamingT()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ip, jp
    integer(kind=4) :: alpha
    
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