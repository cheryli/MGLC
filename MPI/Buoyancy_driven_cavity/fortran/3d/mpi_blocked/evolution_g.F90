subroutine collisionT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !------------------------
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: q(0:6)

    !$omp parallel do default(none) shared(g,g_post,u,v,w,T) private(i,j,k,alpha,n,neq,q,n_post) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
            
    n(0) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(1) = g(1,i,j,k)-g(2,i,j,k)
    n(2) = g(3,i,j,k)-g(4,i,j,k)
    n(3) = g(5,i,j,k)-g(6,i,j,k)
    n(4) = -6.0d0*g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(5) = 2.0d0*g(1,i,j,k)+2.0d0*g(2,i,j,k)-g(3,i,j,k)-g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
    n(6) = g(3,i,j,k)+g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
        
            neq(0) = T(i,j,k)
            neq(1) = T(i,j,k)*u(i,j,k)
            neq(2) = T(i,j,k)*v(i,j,k)
            neq(3) = T(i,j,k)*w(i,j,k)
            neq(4) = T(i,j,k)*paraA
            neq(5) = 0.0d0
            neq(6) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qd
            q(4) = Qnu
            q(5) = Qnu
            q(6) = Qnu
        
            do alpha=0,6
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j,k) = n_post(0)/7.0d0-n_post(4)/7.0d0
    g_post(1,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(2,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(3,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(4,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(5,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
    g_post(6,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
        
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collisionT


subroutine streamingT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(g,g_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,6
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    g(alpha,i,j,k) = g_post(alpha,ip,jp,kp)
                    
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine streamingT


subroutine bouncebackT()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef TopBottomPlatesAdiabatic
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                g(5,i,j,1) = g_post(6,i,j,1)
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                g(6,i,j,nz) = g_post(5,i,j,nz)
            enddo
        enddo
    endif
#endif

#ifdef TopBottomPlatesConstT
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                g(5,i,j,1) = -g_post(6,i,j,1)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                g(6,i,j,nz) = -g_post(5,i,j,nz)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif


#ifdef LeftRightWallsConstT
    ! Left side (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                g(3,i,1,k) = -g_post(4,i,1,k)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                g(4,i,ny,k) = -g_post(3,i,ny,k)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif

#ifdef LeftRightWallsAdiabatic
    ! Left side (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                g(3,i,1,k) = g_post(4,i,1,k)
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                g(4,i,ny,k) = g_post(3,i,ny,k)
            enddo
        enddo
    endif
#endif
    
    
#ifdef BackFrontWallsAdiabatic
    ! Back side (i=1)
    if (coords(0) == 0) then
        do k=1,nz
            do j=1,ny
                g(1,1,j,k) = g_post(2,1,j,k)
            enddo
        enddo
    endif

    !Front side (i=nx)
    if (coords(0) == dims(0) - 1) then
        do k=1,nz
            do j=1,ny          
                g(2,nx,j,k) = g_post(1,nx,j,k)
            enddo
        enddo
    endif

#endif

    return
end subroutine bouncebackT




subroutine macroT()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(g,T) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                T(i,j,k) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macroT